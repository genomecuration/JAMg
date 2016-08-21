static char rcsid[] = "$Id: dynprog_simd.c 189207 2016-05-06 23:16:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_simd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For ceil, log, pow */
#include <ctype.h>		/* For tolower */

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

#include "mem.h"
#include "comp.h"
#include "assert.h"


#ifdef HAVE_AVX2
#define _MM_ADD_EPI8(x,y) _mm256_add_epi8(x,y)
#define _MM_ADDS_EPI8(x,y) _mm256_adds_epi8(x,y)
#define _MM_SUBS_EPI8(x,y) _mm256_subs_epi8(x,y)
#define _MM_CMPGT_EPI8(x,y) _mm256_cmpgt_epi8(x,y)
#define _MM_CMPLT_EPI8(x,y) _mm256_cmpgt_epi8(y,x) /* No _mm256_cmplt commands */
#define _MM_MAX_EPI8(x,y) _mm256_max_epi8(x,y)
#define _MM_MIN_EPI8(x,y) _mm256_min_epi8(x,y)
#define _MM_SET1_EPI8(x) _mm256_set1_epi8(x)

#define _MM_ADD_EPI16(x,y) _mm256_add_epi16(x,y)
#define _MM_ADDS_EPI16(x,y) _mm256_adds_epi16(x,y)
#define _MM_SUBS_EPI16(x,y) _mm256_subs_epi16(x,y)
#define _MM_CMPGT_EPI16(x,y) _mm256_cmpgt_epi16(x,y)
#define _MM_CMPLT_EPI16(x,y) _mm256_cmpgt_epi16(y,x) /* No _mm256_cmplt commands */
#define _MM_MAX_EPI16(x,y) _mm256_max_epi16(x,y)
#define _MM_MIN_EPI16(x,y) _mm256_min_epi16(x,y)
#define _MM_SET1_EPI16(x) _mm256_set1_epi16(x)

#define _MM_SETZERO_SI _mm256_setzero_si256
/* #define _MM_SLLI_SI(x,y) _mm256_slli_si256(x,y) -- 256-bit version works within 128-bit lanes */
/* #define _MM_SRLI_SI(x,y) _mm256_srli_si256(x,y) -- 256-bit version works within 128-bit lanes */
#define _MM_ANDNOT_SI(x,y) _mm256_andnot_si256(x,y)
#define _MM_OR_SI(x,y) _mm256_or_si256(x,y)
#define _MM_AND_SI(x,y) _mm256_and_si256(x,y)

#elif defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
#define _MM_ADD_EPI8(x,y) _mm_add_epi8(x,y)
#define _MM_ADDS_EPI8(x,y) _mm_adds_epi8(x,y)
#define _MM_SUBS_EPI8(x,y) _mm_subs_epi8(x,y)
#define _MM_CMPGT_EPI8(x,y) _mm_cmpgt_epi8(x,y)
#define _MM_CMPLT_EPI8(x,y) _mm_cmplt_epi8(x,y)
#define _MM_MAX_EPI8(x,y) _mm_max_epi8(x,y)
#define _MM_MIN_EPI8(x,y) _mm_min_epi8(x,y)
#define _MM_SET1_EPI8(x) _mm_set1_epi8(x)

#define _MM_ADD_EPI16(x,y) _mm_add_epi16(x,y)
#define _MM_ADDS_EPI16(x,y) _mm_adds_epi16(x,y)
#define _MM_SUBS_EPI16(x,y) _mm_subs_epi16(x,y)
#define _MM_CMPGT_EPI16(x,y) _mm_cmpgt_epi16(x,y)
#define _MM_CMPLT_EPI16(x,y) _mm_cmplt_epi16(x,y)
#define _MM_MAX_EPI16(x,y) _mm_max_epi16(x,y)
#define _MM_MIN_EPI16(x,y) _mm_min_epi16(x,y)
#define _MM_SET1_EPI16(x) _mm_set1_epi16(x)

#define _MM_SETZERO_SI _mm_setzero_si128
/* #define _MM_SLLI_SI(x,y) _mm_slli_si128(x,y) -- 256-bit version works within 128-bit lanes */
/* #define _MM_SRLI_SI(x,y) _mm_srli_si128(x,y) -- 256-bit version works within 128-bit lanes */
#define _MM_ANDNOT_SI(x,y) _mm_andnot_si128(x,y)
#define _MM_OR_SI(x,y) _mm_or_si128(x,y)
#define _MM_AND_SI(x,y) _mm_and_si128(x,y)
#endif



#define LAZY_INDEL 1		/* Don't advance to next coordinate on final indel, since could go over chromosome bounds. */

/* Row 0 and column 0 directions */
/* Was useful in finding a saturation bug, but can fail because of saturation */
#ifdef CHECK1
#define check1(x) x
#else
#define check1(x)
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Fgap */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Compare SIMD with non-SIMD.  Define in dynprog.h */
#ifdef DEBUG_SIMD
#define debug_simd(x) x
#else
#define debug_simd(x)
#endif

#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* Compare AVX2 with SSE42.  Define in dynprog.h */
#ifdef DEBUG_AVX2
#define debug_avx2(x) x
#else
#define debug_avx2(x)
#endif

/* Checking genomic nt in traceback procedures */
#ifdef DEBUG17
#define debug17(x) x
#else
#define debug17(x)
#endif



#include "complement.h"
#define NEG_INFINITY_DISPLAY -99


/************************************************************************
 *   Debugging procedures
 ************************************************************************/

#ifdef DEBUG15
/* For debugging of SIMD procedures*/
#ifdef HAVE_AVX2
static void
print_vector_8 (__m256i x, int r, int c, char *label) {
  __m256i a[1];
  Score8_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm256_store_si256(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15],
	 s[16],s[17],s[18],s[19],s[20],s[21],s[22],s[23],s[24],s[25],s[26],s[27],s[28],s[29],s[30],s[31]);
  return;
}

static void
print_vector_16 (__m256i x, int r, int c, char *label) {
  __m256i a[1];
  Score16_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm256_store_si256(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15]);
  return;
}

#else
static void
print_vector_8 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score8_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	 r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14],s[15]);
  return;
}

static void
print_vector_16 (__m128i x, int r, int c, char *label) {
  __m128i a[1];
  Score16_T *s = a;

  _mm_lfence();			/* Needed to print correct values */
  _mm_store_si128(a,x);
  printf("%d,%d %s: %d %d %d %d %d %d %d %d\n",r,c,label,s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
  return;
}
#endif
#endif


#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD) || defined(DEBUG2)
static void
Matrix8_print (Score8_T **matrix, int rlength, int glength, char *rsequence,
	       char *gsequence, char *gsequencealt,
	       bool revp, int lband, int uband) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d ",j);
  }
  printf("\n");


  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c ",' ');
	} else {
	  printf("  %c ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("  . ");
      } else if (j > i + uband) {
	printf("  . ");
      } else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	printf("%3d ",NEG_INFINITY_DISPLAY);
      } else {
	printf("%3d ",matrix[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Matrix8_print_ud (Score8_T **matrix, int rlength, int glength, char *rsequence,
		  char *gsequence, char *gsequencealt,
		  bool revp, int band, bool upperp) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c ",' ');
	} else {
	  printf("  %c ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    if (upperp == true) {
      for (j = 0; j <= glength; ++j) {
	if (j < i) {
	  printf("  . ");
	} else if (j > i + band) {
	  printf("  . ");
	} else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	  printf("%3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf("%3d ",matrix[j][i]);
	}
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (i < j) {
	  printf("  . ");
	} else if (i > j + band) {
	  printf("  . ");
	} else if (matrix[i][j] < NEG_INFINITY_DISPLAY) {
	  printf("%3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf("%3d ",matrix[i][j]);
	}
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}


static void
Matrix16_print (Score16_T **matrix, int rlength, int glength, char *rsequence,
		char *gsequence, char *gsequencealt,
		bool revp, int lband, int uband) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  if (rlength >= 100) {
    printf("    ");
  } else {
    printf("   ");		/* For i */
  }
  printf("  ");
  if (glength >= 100) {
    for (j = 0; j <= glength; ++j) {
      printf(" %3d ",j);
    }
  } else {
    for (j = 0; j <= glength; ++j) {
      printf(" %2d ",j);
    }
  }
  printf("\n");

  if (gsequence) {
    if (rlength >= 100) {
      printf("    ");
    } else {
      printf("   ");		/* For i */
    }
    printf("  ");
    if (glength >= 100) {
      for (j = 0; j <= glength; ++j) {
	if (j == 0) {
	  printf("     ");
	} else {
	  printf("   %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
	}
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (j == 0) {
	  printf("    ");
	} else {
	  printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
	}
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    if (rlength >= 100) {
      printf("    ");
    } else {
      printf("   ");		/* For i */
    }
    printf("  ");
    if (glength >= 100) {
      for (j = 0; j <= glength; ++j) {
	if (j == 0) {
	  printf("     ");
	} else {
	  g = revp ? gsequence[-j+1] : gsequence[j-1];
	  g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	  if (g == g_alt) {
	    printf("   %c ",' ');
	  } else {
	    printf("   %c ",g_alt);
	  }
	}
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (j == 0) {
	  printf("    ");
	} else {
	  g = revp ? gsequence[-j+1] : gsequence[j-1];
	  g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	  if (g == g_alt) {
	    printf("  %c ",' ');
	  } else {
	    printf("  %c ",g_alt);
	  }
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    if (rlength >= 100) {
      printf("%3d ",i);
    } else {
      printf("%2d ",i);
    }
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    if (glength >= 100) {
      for (j = 0; j <= glength; ++j) {
	if (j < i - lband) {
	  printf("   . ");
	} else if (j > i + uband) {
	  printf("   . ");
	} else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	  printf(" %3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf(" %3d ",matrix[j][i]);
	}
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (j < i - lband) {
	  printf("  . ");
	} else if (j > i + uband) {
	  printf("  . ");
	} else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	  printf("%3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf("%3d ",matrix[j][i]);
	}
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Matrix16_print_ud (Score16_T **matrix, int rlength, int glength, char *rsequence,
		   char *gsequence, char *gsequencealt,
		   bool revp, int band, bool upperp) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	printf("  %c ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("    ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c ",' ');
	} else {
	  printf("  %c ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    if (upperp == true) {
      for (j = 0; j <= glength; ++j) {
	if (j < i) {
	  printf("  . ");
	} else if (j > i + band) {
	  printf("  . ");
	} else if (matrix[j][i] < NEG_INFINITY_DISPLAY) {
	  printf("%3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf("%3d ",matrix[j][i]);
	}
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (i < j) {
	  printf("  . ");
	} else if (i > j + band) {
	  printf("  . ");
	} else if (matrix[i][j] < NEG_INFINITY_DISPLAY) {
	  printf("%3d ",NEG_INFINITY_DISPLAY);
	} else {
	  printf("%3d ",matrix[i][j]);
	}
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif

#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD) || defined(DEBUG2)
static void
Directions8_print (Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		   int rlength, int glength, char *rsequence, char *gsequence, char *gsequencealt,
		   bool revp, int lband, int uband) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d   ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c   ",' ');
	} else {
	  printf("  %c   ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("     ");
      } else if (j > i + uband) {
	printf("     ");
      } else {
	if (directions_Egap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be HORIZ */
	  printf("H");
	}
	printf("|");
	if (directions_nogap[j][i] == DIAG) {
	  printf("D");
	} else if (directions_nogap[j][i] == HORIZ) {
	  printf("H");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
	printf("|");
	if (directions_Fgap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Directions8_print_ud (Direction8_T **directions_nogap, Direction8_T **directions_Egap,
		      int rlength, int glength, char *rsequence, char *gsequence, char *gsequencealt,
		      bool revp, int band, bool upperp) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d   ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c   ",' ');
	} else {
	  printf("  %c   ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    if (upperp == true) {
      for (j = 0; j <= glength; ++j) {
	if (j < i) {
	  printf("     ");
	} else if (j > i + band) {
	  printf("     ");
	} else {
	  if (directions_Egap[j][i] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	  printf("|");
	  if (directions_nogap[j][i] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	  printf("| ");		/* For Fgap */
	}
	printf(" ");
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	if (i < j) {
	  printf("     ");
	} else if (i > j + band) {
	  printf("     ");
	} else {
	  printf(" |");		/* For Fgap */
	  if (directions_nogap[i][j] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	  printf("|");
	  if (directions_Egap[i][j] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	}
	printf(" ");
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}


static void
Directions16_print (Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		    int rlength, int glength, char *rsequence, char *gsequence, char *gsequencealt,
		    bool revp, int lband, int uband) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %3d  ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c   ",' ');
	} else {
	  printf("  %c   ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    for (j = 0; j <= glength; ++j) {
      if (j < i - lband) {
	printf("     ");
      } else if (j > i + uband) {
	printf("     ");
      } else {
	if (directions_Egap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be HORIZ */
	  printf("H");
	}
	printf("|");
	if (directions_nogap[j][i] == DIAG) {
	  printf("D");
	} else if (directions_nogap[j][i] == HORIZ) {
	  printf("H");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
	printf("|");
	if (directions_Fgap[j][i] == DIAG) {
	  printf("D");
	} else {
	  /* Must be VERT */
	  printf("V");
	}
      }
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Directions16_print_ud (Direction16_T **directions_nogap, Direction16_T **directions_Egap,
		       int rlength, int glength, char *rsequence, char *gsequence, char *gsequencealt,
		       bool revp, int band, bool upperp) {
  int i, j;
  char g, g_alt;

  _mm_lfence();

  /* j */
  printf("   ");		/* For i */
  printf("  ");
  for (j = 0; j <= glength; ++j) {
    printf(" %2d   ",j);
  }
  printf("\n");

  if (gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	printf("  %c   ",revp ? gsequence[-j+1] : gsequence[j-1]);
      }
    }
    printf("\n");
  }

  if (gsequencealt != gsequence) {
    printf("   ");		/* For i */
    printf("  ");
    for (j = 0; j <= glength; ++j) {
      if (j == 0) {
	printf("      ");
      } else {
	g = revp ? gsequence[-j+1] : gsequence[j-1];
	g_alt = revp ? gsequencealt[-j+1] : gsequencealt[j-1];
	if (g == g_alt) {
	  printf("  %c   ",' ');
	} else {
	  printf("  %c   ",g_alt);
	}
      }
    }
    printf("\n");
  }

  for (i = 0; i <= rlength; ++i) {
    printf("%2d ",i);
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? rsequence[-i+1] : rsequence[i-1]);
    }
    if (upperp == true) {
      for (j = 0; j <= glength; ++j) {
	if (j < i) {
	  printf("   ");
	} else if (j > i + band) {
	  printf("   ");
	} else {
	  if (directions_Egap[j][i] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	  printf("|");
	  if (directions_nogap[j][i] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	}
	printf("  ");		/* For Fgap */
	printf(" ");
      }
    } else {
      for (j = 0; j <= glength; ++j) {
	printf("  ");		/* For Fgap */
	if (i < j) {
	  printf("   ");
	} else if (i > j + band) {
	  printf("   ");
	} else {
	  if (directions_nogap[i][j] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	  printf("|");
	  if (directions_Egap[i][j] == DIAG) {
	    printf("D");
	  } else {
	    printf("-");
	  }
	}
	printf(" ");
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif



#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
static void
banded_matrix8_compare (Score8_T **matrix1,
#ifdef DEBUG_AVX2
			Score8_T **matrix2,
#else
			Score32_T **matrix2,
#endif
			int rlength, int glength, int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_8 + 30 && matrix2[c][r] <= NEG_INFINITY_8 + 30) {
	/* Okay: both essentially negative infinity */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix8_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
		      revp,lband,uband);
#ifdef DEBUG_AVX2
	Matrix8_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
		      revp,lband,uband);
#elif defined(DEBUG_SIMD)
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
#endif
	abort();
      }
    }
  }

  return;
}

static void
banded_matrix8_compare_upper (Score8_T **matrix1,
#ifdef DEBUG_AVX2
			      Score8_T **matrix2,
#else
			      Score32_T **matrix2,
#endif
			      int rlength, int glength,
			      int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			      bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_8 + 30 && matrix2[c][r] <= NEG_INFINITY_8 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix8_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			 revp,uband,/*upperp*/true);
#ifdef DEBUG_AVX2
	Matrix8_print_ud(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			 revp,uband,/*upperp*/true);
#else
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,/*lband*/0,uband);
#endif
	abort();
      }
    }
  }

  return;
}

static void
banded_matrix8_compare_lower (Score8_T **matrix1,
#ifdef DEBUG_AVX2
			      Score8_T **matrix2,
#else
			      Score32_T **matrix2,
#endif
			      int rlength, int glength,
			      int lband, char *rsequence, char *gsequence, char *gsequence_alt,
			      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			      bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
#ifdef DEBUG_AVX2
      if (matrix1[r][c] <= NEG_INFINITY_8 + 30 && matrix2[r][c] <= NEG_INFINITY_8 + 30) {
	/* Okay */
      } else if (matrix1[r][c] != matrix2[r][c]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[r][c],matrix2[r][c]);

	Matrix8_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			 revp,lband,/*upperp*/false);
	Matrix8_print_ud(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			 revp,lband,/*upperp*/false);
	abort();
      }
#else
      if (matrix1[r][c] <= NEG_INFINITY_8 + 30 && matrix2[c][r] <= NEG_INFINITY_8 + 30) {
	/* Okay */
      } else if (matrix1[r][c] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[r][c],matrix2[c][r]);

	Matrix8_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			 revp,lband,/*upperp*/false);
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,/*uband*/0);
	abort();
      }
#endif
    }
  }

  return;
}


static void
banded_matrix16_compare (Score16_T **matrix1,
#ifdef DEBUG_AVX2
			 Score16_T **matrix2,
#elif defined(DEBUG_SIMD)
			 Score32_T **matrix2,
#endif
			 int rlength, int glength, int lband, int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			 bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_16 + 30 && matrix2[c][r] <= NEG_INFINITY_16 + 30) {
	/* Okay: both essentially negative infinity */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix16_print(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,lband,uband);
#ifdef DEBUG_AVX2
	Matrix16_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,lband,uband);
#elif defined(DEBUG_SIMD)
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,uband);
#endif
	abort();
      }
    }
  }

  return;
}

static void
banded_matrix16_compare_upper (Score16_T **matrix1,
#ifdef DEBUG_AVX2
			       Score16_T **matrix2,
#else
			       Score32_T **matrix2,
#endif
			       int rlength, int glength,
			       int uband, char *rsequence, char *gsequence, char *gsequence_alt,
			       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			       bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix1[c][r] <= NEG_INFINITY_16 + 30 && matrix2[c][r] <= NEG_INFINITY_16 + 30) {
	/* Okay */
      } else if (matrix1[c][r] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[c][r],matrix2[c][r]);

	Matrix16_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			    revp,uband,/*upperp*/true);
#ifdef DEBUG_AVX2
	Matrix16_print_ud(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			    revp,uband,/*upperp*/true);
#else
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
				 goffset,chroffset,chrhigh,watsonp,revp,/*lband*/0,uband);
#endif
	abort();
      }
    }
  }

  return;
}

static void
banded_matrix16_compare_lower (Score16_T **matrix1,
#ifdef DEBUG_AVX2
			       Score16_T **matrix2,
#else
			       Score32_T **matrix2,
#endif
			       int rlength, int glength,
			       int lband, char *rsequence, char *gsequence, char *gsequence_alt,
			       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			       bool revp) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
#ifdef DEBUG_AVX2
      if (matrix1[r][c] <= NEG_INFINITY_16 + 30 && matrix2[r][c] <= NEG_INFINITY_16 + 30) {
	/* Okay */
      } else if (matrix1[r][c] != matrix2[r][c]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[r][c],matrix2[r][c]);

	Matrix16_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			  revp,lband,/*upperp*/false);
	Matrix16_print_ud(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			  revp,lband,/*upperp*/false);
	abort();
      }
#else
      if (matrix1[r][c] <= NEG_INFINITY_16 + 30 && matrix2[c][r] <= NEG_INFINITY_16 + 30) {
	/* Okay */
      } else if (matrix1[r][c] != matrix2[c][r]) {
	printf("At %d,%d, value %d != value %d\n",r,c,matrix1[r][c],matrix2[c][r]);

	Matrix16_print_ud(matrix1,rlength,glength,rsequence,gsequence,gsequence_alt,
			  revp,lband,/*upperp*/false);
	Dynprog_Matrix32_print(matrix2,rlength,glength,rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp,lband,/*uband*/0);
	abort();
      }
#endif
    }
  }

  return;
}

#endif

#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
static void
banded_directions8_compare_nogap (Score8_T **matrix, Direction8_T **directions1,
#ifdef DEBUG_AVX2
				  Direction8_T **directions2,
#elif defined(DEBUG_SIMD)
				  Direction32_T **directions2,
#endif
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions8_compare_nogap_upper (Score8_T **matrix, Direction8_T **directions1,
#ifdef DEBUG_AVX2
					Direction8_T **directions2,
#elif defined(DEBUG_SIMD)
					Direction32_T **directions2,
#endif
					int rlength, int glength, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions8_compare_nogap_lower (Score8_T **matrix, Direction8_T **directions1,
#ifdef DEBUG_AVX2
					Direction8_T **directions2,
#elif defined(DEBUG_SIMD)
					Direction32_T **directions2,
#endif
					int rlength, int glength, int lband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

#ifdef DEBUG_AVX2
    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[r][c] == 0) {
	if (directions2[r][c] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[r][c] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else {
	if (directions2[r][c] == 0 || directions2[r][c] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}
      }
    }

#else
    for (r = rlo; r <= rhigh; r++) {
      if (matrix[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[r][c] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}
      }
    }
#endif

  }

  return;
}


static void
banded_directions16_compare_nogap (Direction16_T **directions1,
#ifdef DEBUG_AVX2
				   Direction16_T **directions2,
#elif defined(DEBUG_SIMD)
				   Direction32_T **directions2,
#endif
				   int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_nogap_upper (Direction16_T **directions1,
#ifdef DEBUG_AVX2
					 Direction16_T **directions2,
#else
					 Direction32_T **directions2,
#endif
					 int rlength, int glength, int uband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_nogap_lower (Direction16_T **directions1,
#ifdef DEBUG_AVX2
					 Direction16_T **directions2,
#else
					 Direction32_T **directions2,
#endif
					 int rlength, int glength, int lband) {
  int r, c, rlo, rhigh;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = rlo; r <= rhigh; r++) {
#ifdef DEBUG_AVX2
      if (directions1[r][c] == 0) {
	if (directions2[r][c] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[r][c] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else {
	if (directions2[r][c] == 0 || directions2[r][c] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}
      }
#else
      if (directions1[r][c] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, nogap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}
      }
#endif
    }
  }

  return;
}
#endif

#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
static void
banded_directions8_compare_Egap (Score8_T **matrix1, Direction8_T **directions1,
#ifdef DEBUG_AVX2
				 Direction8_T **directions2,
#else
				 Direction32_T **directions2,
#endif
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions8_compare_Egap_upper (Score8_T **matrix1, Direction8_T **directions1,
#ifdef DEBUG_AVX2
				       Direction8_T **directions2,
#else
				       Direction32_T **directions2,
#endif
				       int rlength, int glength, int uband) {
  int r, c, rlo, rhigh, last_check;

  return;
  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions8_compare_Egap_lower (Score8_T **matrix1, Direction8_T **directions1,
#ifdef DEBUG_AVX2
				       Direction8_T **directions2,
#else
				       Direction32_T **directions2,
#endif
				       int rlength, int glength, int lband) {
  int r, c, rlo, rhigh, last_check;

  return;
  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
#ifdef DEBUG_AVX2
      if (matrix1[r][c] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[r][c] == 0) {
	if (directions2[r][c] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[r][c] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else {
	if (directions2[r][c] == 0 || directions2[r][c] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}
      }
#else
      if (matrix1[r][c] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[r][c] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else if (directions1[r][c] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}
      }
#endif
    }
  }

  return;
}


static void
banded_directions16_compare_Egap (Direction16_T **directions1,
#ifdef DEBUG_AVX2
				  Direction16_T **directions2,
#else
				  Direction32_T **directions2,
#endif
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, last_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Egap_upper (Direction16_T **directions1,
#ifdef DEBUG_AVX2
					Direction16_T **directions2,
#else
					Direction32_T **directions2,
#endif
					int rlength, int glength, int uband) {
  int r, c, rlo, rhigh, last_check;

  return;
  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      rlo = 1;
    };

    if ((rhigh = c) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Egap_lower (Direction16_T **directions1,
#ifdef DEBUG_AVX2
					Direction16_T **directions2,
#else
					Direction32_T **directions2,
#endif
					int rlength, int glength, int lband) {
  int r, c, rlo, rhigh, last_check;

  return;
  for (c = 1; c <= glength; c++) {
    if ((rlo = c) < 1) {
      rlo = 1;
    };

    if ((rhigh = c + lband) <= rlength) {
      /* Don't check rhigh.  Egap direction derives from a comparison
	 of NEG_INFINITY values, and we should never reach here from
	 directions_nogap anyway. */
      last_check = rhigh - 1;

    } else {
      /* Do check rhigh, which contains instructions for the bottom row */
      rhigh = rlength;
      last_check = rhigh;
    }

    for (r = rlo; r <= last_check; r++) {
#ifdef DEBUG_AVX2
      if (directions1[r][c] == 0) {
	if (directions2[r][c] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}
      } else if (directions1[r][c] == 1) {
	if (directions2[r][c] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}

      } else {
	if (directions2[r][c] == 0 || directions2[r][c] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[r][c]);
	  abort();
	}
      }
#else
      if (directions1[r][c] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}
      } else if (directions1[r][c] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Egap dir %d != dir %d\n",r,c,directions1[r][c],directions2[c][r]);
	  abort();
	}
      }
#endif
    }
  }

  return;
}
#endif


#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
static void
banded_directions8_compare_Fgap (Score8_T **matrix1, Direction8_T **directions1,
#ifdef DEBUG_AVX2
				 Direction8_T **directions2,
#else
				 Direction32_T **directions2,
#endif
				 int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (matrix1[c][r] < NEG_INFINITY_8 + 30) {
	/* Don't check */

      } else if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}

      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d.  Score is %d\n",
		 r,c,directions1[c][r],directions2[c][r],matrix1[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}

static void
banded_directions16_compare_Fgap (Direction16_T **directions1,
#ifdef DEBUG_AVX2
				  Direction16_T **directions2,
#else
				  Direction32_T **directions2,
#endif
				  int rlength, int glength, int lband, int uband) {
  int r, c, rlo, rhigh, first_check;

  for (c = 1; c <= glength; c++) {
    if ((rlo = c - uband) < 1) {
      first_check = rlo = 1;
    } else {
      first_check = rlo + 1;
    }

    if ((rhigh = c + lband) > rlength) {
      rhigh = rlength;
    }

    for (r = first_check; r <= rhigh; r++) {
      if (directions1[c][r] == 0) {
	if (directions2[c][r] == 0) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      } else if (directions1[c][r] == 1) {
	if (directions2[c][r] == 1) {
	} else {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}

      } else {
	if (directions2[c][r] == 0 || directions2[c][r] == 0) {
	  printf("At %d,%d, Fgap dir %d != dir %d\n",r,c,directions1[c][r],directions2[c][r]);
	  abort();
	}
      }
    }
  }

  return;
}
#endif


/************************************************************************
 *   End of debugging procedures
 ************************************************************************/



#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score8_T **
aligned_score8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = ptr;	   /* Want aligned row to be r = 0, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	   /* Want aligned row to be r = 0, 16, ... */
  }
#if defined(DEBUG2) && (defined(DEBUG_AVX2) || defined(DEBUG_SIMD))
  memset((void *) matrix[0],0,(glength+1)*rlength*sizeof(Score8_T));
#endif

  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score8_T **
aligned_directions8_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = ptr;	   /* Want aligned row to be r = 0, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	   /* Want aligned row to be r = 0, 16, ... */
  }
#if defined(DEBUG2) && (defined(DEBUG_AVX2) || defined(DEBUG_SIMD))
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*rlength*sizeof(Score8_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score8_T **
aligned_directions8_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score8_T **matrix, *ptr;
  int c;

  matrix = (Score8_T **) ptrs;

  ptr = (Score8_T *) space;
  matrix[0] = ptr;	/* Want aligned row to be r = 0, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	/* Want aligned row to be r = 0, 16, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*rlength*sizeof(Score8_T));

  return matrix;
}



/* Makes a matrix of dimensions 0..rlength x 0..glength inclusive */
static Score16_T **
aligned_score16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],0,(glength+1)*rlength*sizeof(Score16_T));
#endif

  return matrix;
}

/* No initialization to DIAG (0), for directions_Egap and directions_nogap */
static Score16_T **
aligned_directions16_alloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  }
#ifdef DEBUG2
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*rlength*sizeof(Score16_T));
#endif

  return matrix;
}

/* Initialization to DIAG (0), for directions_Fgap */
static Score16_T **
aligned_directions16_calloc (int rlength, int glength, void **ptrs, void *space) {
  Score16_T **matrix, *ptr;
  int c;

  matrix = (Score16_T **) ptrs;

  ptr = (Score16_T *) space;
  matrix[0] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  for (c = 1; c <= glength; c++) {
    ptr += rlength;
    matrix[c] = ptr;	/* Want aligned row to be r = 0, 8, 16, ... */
  }
  memset((void *) matrix[0],/*DIAG*/0,(glength+1)*rlength*sizeof(Score16_T));

  return matrix;
}
#endif


#define T Dynprog_T


#ifdef DEBUG_AVX2
Score8_T **
Dynprog_simd_8_nonavx2 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
			Direction8_T ***directions_Fgap,
			T this, char *rsequence, char *gsequence, char *gsequence_alt,
			int rlength, int glength,
			int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			Mismatchtype_T mismatchtype, int open, int extend,
			int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF;	/* Need to have the ability to go past NEG_INFINITY */
  Score8_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
#ifndef HAVE_SSE4_1
  __m128i pairscores_best, all_128;
#endif
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, extend_ladder, extend_chunk, complement_dummy;
  __m128i dir_horiz;
  __m128i ramp, ramp_chunk, lband_vector, filter, ones;
  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score8_T *pairscores_col0;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr, pairscore, pairscore0;
  Pairdistance_T **pairdistance_array_type;


  debug2(printf("Dynprog_simd_8.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8.  jump_late_p %d, open %d, extend %d, lband %d, uband %d\n",
		 jump_late_p,open,extend,lband,uband));

  rlength_ceil = (int) ((rlength + SIMD_NCHARS_NONAVX2)/SIMD_NCHARS_NONAVX2) * SIMD_NCHARS_NONAVX2;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Need to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("Dynprog_simd_8: "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned_std.one.matrix_ptrs,this->aligned_std.one.matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned_std.one.directions_ptrs_0,this->aligned_std.one.directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned_std.one.directions_ptrs_1,this->aligned_std.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions8_calloc(rlength_ceil,glength,
						this->aligned_std.one.directions_ptrs_2,this->aligned_std.one.directions_space_2);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NCHARS_NONAVX2 && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);

  lband_ceil = (int) ((lband + SIMD_NCHARS_NONAVX2)/SIMD_NCHARS_NONAVX2) * SIMD_NCHARS_NONAVX2;
  pairscores_col0 = (Score8_T *) _mm_malloc(lband_ceil * sizeof(Score8_T),16);


#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
#ifdef HAVE_SSE4_1
  pairscores_col0[0] = (Score8_T) 0;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) 0;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) NEG_INFINITY_8;
  }
#endif
#else
  pairscores_col0[0] = (Score8_T) 0+128;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) 0+128;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) NEG_INFINITY_8+128;
  }
#endif
#endif


  /* Row 0 */
  r = 0; na1 = 'N';
  pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
#endif

  complement_dummy = _mm_set1_epi8(-1);

  FF = (int *) MALLOCA((glength + 1) * sizeof(int));

  gap_open = _mm_set1_epi8((Score8_T) open);
  gap_extend = _mm_set1_epi8((Score8_T) extend);

#ifndef INFINITE_INITIAL_GAP_PENALTY
  lband_vector = _mm_set1_epi8(lband);
  ramp = _mm_setr_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
  extend_ladder = _mm_setr_epi8(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend,
				9*extend,10*extend,11*extend,12*extend,13*extend,14*extend,15*extend,16*extend);
  ramp_chunk = _mm_set1_epi8(SIMD_NCHARS_NONAVX2);
  extend_chunk = _mm_set1_epi8(SIMD_NCHARS_NONAVX2*extend);
#endif

  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NCHARS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8-open);
	filter = _mm_cmpgt_epi8(ramp,lband_vector);
	H_nogap_r = _mm_and_si128(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _mm_adds_epi8(E_r_gap,gap_open);

	ramp = _mm_adds_epi8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi8(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E >= H.  To fill in first column of each
	   row block with non-diags, make E == H. */
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8);
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   equal to E. */
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* To compensate for T1 = H + open */
	filter = _mm_cmpgt_epi8(ramp,lband_vector);
	H_nogap_r = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),_mm_and_si128(filter,H_nogap_r));
	E_r_gap = _mm_adds_epi8(H_nogap_r,gap_open);
	ramp = _mm_adds_epi8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi8(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
      } else {
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8);
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_set1_epi8(0);
	  } else {
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _mm_set1_epi8(0);
#else
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	  }
	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi8(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
#ifndef HAVE_SSE4_1
	    pairscore -= 128;
#endif
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_8) {
	      score_column[rhigh_calc] = NEG_INFINITY_8; /* Saturation */
	    } else if (score > POS_INFINITY_8) {
	      /* Should never get here, because we limit size of matrix using 8-bit quantities */
	      score_column[rhigh_calc] = POS_INFINITY_8; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score8_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }

  } else {
    /* jump early */
#if defined(ZERO_INITIAL_GAP_PENALTY) || defined(INFINITE_INITIAL_GAP_PENALTY)
    /* No need for ones */
#else
    ones = _mm_set1_epi8(1);
#endif
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NCHARS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8-open);
	filter = _mm_cmpgt_epi8(ramp,lband_vector);
	H_nogap_r = _mm_and_si128(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _mm_adds_epi8(E_r_gap,gap_open);

	ramp = _mm_adds_epi8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi8(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E > H.  To fill in first column of each
	   row block with non-diags, make E > H. */
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8+1);
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   less than E. */
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open+1); /* To compensate for T1 = H + open */
	filter = _mm_cmpgt_epi8(ramp,lband_vector);
	H_nogap_r = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),_mm_and_si128(filter,H_nogap_r));
	E_r_gap = _mm_adds_epi8(H_nogap_r,gap_open);
	H_nogap_r = _mm_subs_epi8(H_nogap_r,ones); /* To ensure H < E */
	ramp = _mm_adds_epi8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi8(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
      } else {
	E_r_gap = _mm_set1_epi8(NEG_INFINITY_8+1);
	H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_set1_epi8(0);
	  } else {
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _mm_set1_epi8(0);
#else
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	  }
	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi8(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
#else
	/* Compare H + open with vert */
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi8(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
#ifndef HAVE_SSE4_1
	    pairscore -= 128;
#endif
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_8) {
	      score_column[rhigh_calc] = NEG_INFINITY_8; /* Saturation */
	    } else if (score > POS_INFINITY_8) {
	      /* Should never get here, because we limit size of matrix using 8-bit quantities */
	      score_column[rhigh_calc] = POS_INFINITY_8; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score8_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	  debug3(printf("LAST_NOGAP gets score_column[%d-1], or %d\n",rlo_calc,last_nogap));
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    debug3(printf("Stored at score_column[%d]: %d\n",r,(Score8_T) score_column[r]));
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
  (*directions_Fgap)[0][1] = VERT;
#endif  

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8\n");
  Matrix8_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		revp,lband,uband);
  Directions8_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			    rlength,glength,rsequence,gsequence,gsequence_alt,
			    revp,lband,uband);
#endif
  
#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Fgap)[0][r] != DIAG);
    assert((*directions_nogap)[0][r] != DIAG);
  }
#endif

  FREEA(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif



#if defined(HAVE_SSE2)
/* Modified from Dynprog_simd_8_upper.  Operates by columns. */
Score8_T **
Dynprog_simd_8 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
		Direction8_T ***directions_Fgap,
		T this, char *rsequence, char *gsequence, char *gsequence_alt,
		int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		Mismatchtype_T mismatchtype, int open, int extend,
		int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF;	/* Need to have the ability to go past NEG_INFINITY */
  Score8_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i pairscores_std, pairscores_alt;
  __m256i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_horiz;
#if defined(ZERO_INITIAL_GAP_PENALTY)
  __m256i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter;
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
#else
  __m256i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter, ones;
#endif

#else
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
#if defined(ZERO_INITIAL_GAP_PENALTY)
  __m128i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter;
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
#else
  __m128i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter, ones;
#endif
#endif
#ifndef HAVE_SSE4_1
  __m128i pairscores_best, all_128;
#endif
  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score8_T *pairscores_col0;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr, pairscore, pairscore0;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG_AVX2
  Score8_T **matrix_std;
  Direction8_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#endif


  debug2(printf("Dynprog_simd_8.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8.  jump_late_p %d, open %d, extend %d, lband %d, uband %d\n",
		 jump_late_p,open,extend,lband,uband));

  rlength_ceil = (int) ((rlength + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("Dynprog_simd_8: "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned.one.matrix_ptrs,this->aligned.one.matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_0,this->aligned.one.directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned.one.directions_ptrs_1,this->aligned.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions8_calloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_2,this->aligned.one.directions_space_2);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NCHARS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);

  lband_ceil = (int) ((lband + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;
  pairscores_col0 = (Score8_T *) _mm_malloc(lband_ceil * sizeof(Score8_T),ALIGN_SIZE);


#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
#ifdef HAVE_SSE4_1
  pairscores_col0[0] = (Score8_T) 0;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) 0;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) NEG_INFINITY_8;
  }
#endif
#else
  pairscores_col0[0] = (Score8_T) 0+128;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) 0+128;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score8_T) NEG_INFINITY_8+128;
  }
#endif
#endif


  /* Row 0 */
  r = 0; na1 = 'N';
  pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
#endif

  complement_dummy = _MM_SET1_EPI8(-1);

  FF = (int *) MALLOCA((glength + 1) * sizeof(int));

  gap_open = _MM_SET1_EPI8((Score8_T) open);
  gap_extend = _MM_SET1_EPI8((Score8_T) extend);


#ifndef INFINITE_INITIAL_GAP_PENALTY
#ifdef HAVE_AVX2
  ramp = _mm256_setr_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
  extend_ladder = _mm256_setr_epi8(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend,
				   9*extend,10*extend,11*extend,12*extend,13*extend,14*extend,15*extend,16*extend);
#else
  ramp = _mm_setr_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
  extend_ladder = _mm_setr_epi8(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend,
				9*extend,10*extend,11*extend,12*extend,13*extend,14*extend,15*extend,16*extend);
#endif
  lband_vector = _MM_SET1_EPI8(lband);
  ramp_chunk = _MM_SET1_EPI8(SIMD_NCHARS);
  extend_chunk = _MM_SET1_EPI8(SIMD_NCHARS*extend);
#endif

  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8-open);
	filter = _MM_CMPGT_EPI8(ramp,lband_vector);
	H_nogap_r = _MM_AND_SI(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _MM_OR_SI(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _MM_ADDS_EPI8(E_r_gap,gap_open);

	ramp = _MM_ADDS_EPI8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI8(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E >= H.  To fill in first column of each
	   row block with non-diags, make E == H. */
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8);
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   equal to E. */
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* To compensate for T1 = H + open */
	filter = _MM_CMPGT_EPI8(ramp,lband_vector);
	H_nogap_r = _MM_OR_SI(_MM_ANDNOT_SI(filter,extend_ladder),_MM_AND_SI(filter,H_nogap_r));
	E_r_gap = _MM_ADDS_EPI8(H_nogap_r,gap_open);
	ramp = _MM_ADDS_EPI8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI8(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
      } else {
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8);
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  if (rlo == 0) {
	    X_prev_nogap = _mm256_setzero_si256();
	  } else {
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
	  }
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_setzero_si128();
	  } else {
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

#ifdef HAVE_AVX2
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#endif
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_CHAR_INSERT);
	  }

#else
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif	    
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
	  }
#endif
	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPLT_EPI8(E_r_gap,T1); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MAX_EPI8(E_r_gap, T1); /* Compare H + open with vert */
#else
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _MM_ADDS_EPI8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_r,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_r,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi8(H_nogap_r, _mm256_max_epi8(pairscores_std,pairscores_alt));
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPLT_EPI8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _MM_MAX_EPI8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif

	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
#ifndef HAVE_SSE4_1
	    pairscore -= 128;
#endif
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_8) {
	      score_column[rhigh_calc] = NEG_INFINITY_8; /* Saturation */
	    } else if (score > POS_INFINITY_8) {
	      /* Should never get here, because we limit size of matrix using 8-bit quantities */
	      score_column[rhigh_calc] = POS_INFINITY_8; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score8_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
#ifdef HAVE_AVX2
	H_nogap_r = _mm256_load_si256((__m256i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#else
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#endif
      }
    }

  } else {
    /* jump early */
#if defined(ZERO_INITIAL_GAP_PENALTY) || defined(INFINITE_INITIAL_GAP_PENALTY)
    /* No need for ones */
#else
    ones = _MM_SET1_EPI8(1);
#endif
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8-open);
	filter = _MM_CMPGT_EPI8(ramp,lband_vector);
	H_nogap_r = _MM_AND_SI(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _MM_OR_SI(_MM_ANDNOT_SI(filter,extend_ladder),H_nogap_r);
	E_r_gap = _MM_ADDS_EPI8(E_r_gap,gap_open);

	ramp = _MM_ADDS_EPI8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI8(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E > H.  To fill in first column of each
	   row block with non-diags, make E > H. */
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8+1);
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   less than E. */
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open+1); /* To compensate for T1 = H + open */
	filter = _MM_CMPGT_EPI8(ramp,lband_vector);
	H_nogap_r = _MM_OR_SI(_MM_ANDNOT_SI(filter,extend_ladder),_MM_AND_SI(filter,H_nogap_r));
	E_r_gap = _MM_ADDS_EPI8(H_nogap_r,gap_open);
	H_nogap_r = _MM_SUBS_EPI8(H_nogap_r,ones); /* To ensure H < E */
	ramp = _MM_ADDS_EPI8(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI8(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
      } else {
	E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8+1);
	H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  if (rlo == 0) {
	    X_prev_nogap = _mm256_setzero_si256();
	  } else {
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
	  }
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_setzero_si128();
	  } else {
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

#ifdef HAVE_AVX2
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#endif
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_CHAR_INSERT);
	  }

#else
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	  } else {
	    /* second or greater block of 8 */
	    X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
	  }
#endif

	}

	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPGT_EPI8(E_r_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MAX_EPI8(E_r_gap, T1); /* Compare H + open with vert */
#else
	/* Compare H + open with vert */
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
#endif
	E_r_gap = _MM_ADDS_EPI8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_r,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_r,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi8(H_nogap_r, _mm256_max_epi8(pairscores_std,pairscores_alt));
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPGT_EPI8(E_r_gap,H_nogap_r); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _MM_MAX_EPI8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
#ifndef HAVE_SSE4_1
	    pairscore -= 128;
#endif
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_8) {
	      score_column[rhigh_calc] = NEG_INFINITY_8; /* Saturation */
	    } else if (score > POS_INFINITY_8) {
	      /* Should never get here, because we limit size of matrix using 8-bit quantities */
	      score_column[rhigh_calc] = POS_INFINITY_8; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score8_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	  debug3(printf("LAST_NOGAP gets score_column[%d-1], or %d\n",rlo_calc,last_nogap));
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_8) ? NEG_INFINITY_8 : (Score8_T) c_gap; /* Saturation */
	    debug3(printf("Stored at score_column[%d]: %d\n",r,(Score8_T) score_column[r]));
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
#ifdef HAVE_AVX2
	H_nogap_r = _mm256_load_si256((__m256i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#else
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#endif
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
  (*directions_Fgap)[0][1] = VERT;
#endif  

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8\n");
  Matrix8_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		revp,lband,uband);
  Directions8_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			    rlength,glength,rsequence,gsequence,gsequence_alt,
			    revp,lband,uband);
#endif
  
#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Fgap)[0][r] != DIAG);
    assert((*directions_nogap)[0][r] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_8_nonavx2(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				      this,rsequence,gsequence,gsequence_alt,
				      rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				      open,extend,lband,uband,jump_late_p,revp);
#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,lband,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_8,
				/*upperp*/true,/*lowerp*/true);
#endif

#ifdef DEBUG2
  printf("Banded %s\n",revp ? "rev" : "fwd");
  Matrix8_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,revp,lband,uband);
  Directions8_print(*directions_nogap,*directions_Egap,*directions_Fgap,
		    rlength,glength,rsequence,gsequence,gsequence_alt,revp,lband,uband);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix8_compare(matrix,matrix_std,rlength,glength,lband,uband,
			 rsequence,gsequence,gsequence_alt,
			 goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions8_compare_nogap(matrix,*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);
  banded_directions8_compare_Egap(matrix,*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);
  banded_directions8_compare_Fgap(matrix,*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  FREEA(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef DEBUG_AVX2
/* Designed for computation above the main diagonal, so no F loop or bottom masking needed */
/* Operates by columns */
Score8_T **
Dynprog_simd_8_upper_nonavx2 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
			      T this, char *rsequence, char *gsequence, char *gsequence_alt,
			      int rlength, int glength,
			      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			      Mismatchtype_T mismatchtype, int open, int extend,
			      int uband, bool jump_late_p, bool revp) {
  Score8_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
#ifdef HAVE_SSE4_1
  __m128i E_infinity;
#else
  __m128i E_infinity_plus_128;
  __m128i pairscores_best, all_128;
#endif
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
  int rlength_ceil, r, c;
  int rlo, rhigh;
  int na1, na2, na2_alt;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;


  debug2(printf("Dynprog_simd_8_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  rlength_ceil = (int) ((rlength + SIMD_NCHARS_NONAVX2)/SIMD_NCHARS_NONAVX2) * SIMD_NCHARS_NONAVX2;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("compute_scores_simd_8_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned_std.two.upper_matrix_ptrs,this->aligned_std.two.upper_matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned_std.two.upper_directions_ptrs_0,this->aligned_std.two.upper_directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned_std.two.upper_directions_ptrs_1,this->aligned_std.two.upper_directions_space_1);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NCHARS_NONAVX2 && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
  r = 0; na1 = 'N';
  pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
#endif

  complement_dummy = _mm_set1_epi8(-1);

  gap_open = _mm_set1_epi8((Score8_T) open);
  gap_extend = _mm_set1_epi8((Score8_T) extend);


#ifdef HAVE_SSE4_1
  E_infinity = _mm_set1_epi8(POS_INFINITY_8);
#else
  E_infinity_plus_128 = _mm_set1_epi8(POS_INFINITY_8+128);
#endif
  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NCHARS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H .  To fill in first column of each
	 row block with non-diags, could make E == H.  But irrelevant,
	 because these are below the diagonal. */
      E_mask = _mm_set1_epi8(1);

      E_r_gap = _mm_set1_epi8(NEG_INFINITY_8);
      H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];
	
	if (c == 0) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	}

	debug15(print_vector_8(E_mask,rlo,c,"E_mask"));
#ifdef HAVE_SSE4_1
	E_r_gap = _mm_min_epi8(E_r_gap,_mm_add_epi8(E_mask,E_infinity));
#else
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi8(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_min_epi8(E_r_gap,_mm_add_epi8(E_mask,E_infinity));
#elif 1
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_r_gap = _mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128));
	E_r_gap = _mm_add_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(E_r_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	debug15(print_vector_8(pairscores_best,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));

#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
      }
    }

  } else {
    /* jump early */
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NCHARS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with non-diags, could make E > H.  But irrelevant,
	 because these are below the diagonal. */
      E_mask = _mm_set1_epi8(1);

      E_r_gap = _mm_set1_epi8(NEG_INFINITY_8+1);
      H_nogap_r = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */;
	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	}

	debug15(print_vector_8(E_mask,rlo,c,"E_mask"));
#ifdef HAVE_SSE4_1
	E_r_gap = _mm_min_epi8(E_r_gap,_mm_add_epi8(E_mask,E_infinity));
#else
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi8(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

	/* Compare H + open with vert */
#ifdef HAVE_SSE4_1
	E_r_gap = _mm_max_epi8(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_min_epi8(E_r_gap,_mm_add_epi8(E_mask,E_infinity));
#elif 1
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_r_gap = _mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128));
	E_r_gap = _mm_add_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(E_r_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	debug15(print_vector_8(pairscores_best,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi8(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _mm_max_epi8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);

	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
      }
    }
  }

#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8_upper\n");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		   revp,uband,/*upperp*/true);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,uband,/*upperp*/true);
#endif
  
#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


/* E_mask works at the wraparound from POS_INFINITY to NEG_INFINITY.
   It is designed to prevent a horizontal/vertical jump into the empty
   triangle, by setting horizontal/vertical scores to be as small as
   possible, e.g., -128.  However, it is possible that H is also -128,
   so we still need to fix the directions along the main diagonal.

   E_mask shifted:    0    0    0    0    1    1    1    1
   add E_infinity:  127  127  127  127  127  127  127  127
   resulting mask:  127  127  127  127 -128 -128 -128 -128

   To deal with non-SSE4.1 systems, which lack _mm_min_epi8, we need
   too add 128 to E and mask, then take _mm_min_epu8, then subtract
   128, as follows:

   E_mask shifted:    0    0    0    0    1    1    1    1
   add E_inf+128:   255  255  255  255  255  255  255  255
   resulting mask:  255  255  255  255    0    0    0    0
   (compare w/E+128)

*/


#ifdef HAVE_SSE2
/* Designed for computation above the main diagonal, so no F loop or bottom masking needed */
/* Operates by columns */
Score8_T **
Dynprog_simd_8_upper (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
		      T this, char *rsequence, char *gsequence, char *gsequence_alt,
		      int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		      Mismatchtype_T mismatchtype, int open, int extend,
		      int uband, bool jump_late_p, bool revp) {
  Score8_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i E_infinity;
#elif defined(HAVE_SSE4_1)
  __m128i E_infinity;
#else
  __m128i E_infinity_plus_128;
  __m128i pairscores_best, all_128;
#endif
#ifdef HAVE_AVX2
  __m256i pairscores_std, pairscores_alt;
  __m256i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_horiz;
  Score8_T save;
#else
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
#endif
  int rlength_ceil, r, c;
  int rlo, rhigh;
  int na1, na2, na2_alt;
  Score8_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG_AVX2
  Score8_T **matrix_std;
  Direction8_T **directions_nogap_std, **directions_Egap_std;
  char na2_single;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif


  debug2(printf("Dynprog_simd_8_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  rlength_ceil = (int) ((rlength + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("compute_scores_simd_8_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score8_alloc(rlength_ceil,glength,
				this->aligned.two.upper_matrix_ptrs,this->aligned.two.upper_matrix_space);
  *directions_nogap = aligned_directions8_alloc(rlength_ceil,glength,
						this->aligned.two.upper_directions_ptrs_0,this->aligned.two.upper_directions_space_0);
  *directions_Egap = aligned_directions8_alloc(rlength_ceil,glength,
					       this->aligned.two.upper_directions_ptrs_1,this->aligned.two.upper_directions_space_1);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NCHARS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[1] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[2] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[3] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[4] = (Score8_T *) _mm_malloc(rlength_ceil * sizeof(Score8_T),ALIGN_SIZE);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
  r = 0; na1 = 'N';
  pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score8_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score8_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score8_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score8_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score8_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score8_T));
#endif

  complement_dummy = _MM_SET1_EPI8(-1);

  gap_open = _MM_SET1_EPI8((Score8_T) open);
  gap_extend = _MM_SET1_EPI8((Score8_T) extend);

#ifdef HAVE_SSE4_1
  E_infinity = _MM_SET1_EPI8(POS_INFINITY_8);
#else
  E_infinity_plus_128 = _mm_set1_epi8(POS_INFINITY_8+128);
#endif
  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H .  To fill in first column of each
	 row block with non-diags, could make E == H.  But irrelevant,
	 because these are below the diagonal. */
      E_mask = _MM_SET1_EPI8(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8);
      H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];
	
	if (c == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	} else {
	  /* second or greater block of 8 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	}

	debug15(print_vector_8(E_mask,rlo,c,"E_mask"));
#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MIN_EPI8(E_r_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#else
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPLT_EPI8(E_r_gap,T1); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MAX_EPI8(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _MM_MIN_EPI8(E_r_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#elif 1
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_r_gap = _mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128));
	E_r_gap = _mm_add_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(E_r_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_r,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_r,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _MM_ADDS_EPI8(H_nogap_r, _MM_MAX_EPI8(pairscores_std,pairscores_alt));
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	debug15(print_vector_8(pairscores_best,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPLT_EPI8(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));

#ifdef HAVE_SSE4_1
	H_nogap_r = _MM_MAX_EPI8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif


	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi8(E_mask,15);
	E_mask = _mm256_slli_si256(E_mask,ONE_CHAR);
	E_mask = _mm256_insert_epi8(E_mask,save,16);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
#endif
      }
    }

  } else {
    /* jump early */
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NCHARS) {
      if ((rhigh = rlo + SIMD_NCHARS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with non-diags, could make E > H.  But irrelevant,
	 because these are below the diagonal. */
      E_mask = _MM_SET1_EPI8(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_r_gap = _MM_SET1_EPI8(NEG_INFINITY_8+1);
      H_nogap_r = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */;
	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	} else {
	  /* second or greater block of 8 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	}

	debug15(print_vector_8(E_mask,rlo,c,"E_mask"));
#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MIN_EPI8(E_r_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#else
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPGT_EPI8(E_r_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(T1,rlo,c,"T1"));
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_Egap"));

	/* Compare H + open with vert */
#ifdef HAVE_SSE4_1
	E_r_gap = _MM_MAX_EPI8(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _MM_MIN_EPI8(E_r_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#elif 1
	E_r_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_r_gap = _mm_adds_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_r_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_r_gap = _mm_max_epu8(_mm_add_epi8(E_r_gap, all_128), _mm_add_epi8(T1, all_128));
	E_r_gap = _mm_add_epi8(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_sub_epi8(_mm_min_epu8(E_r_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_r,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_r,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_CHAR_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_CHAR);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _MM_ADDS_EPI8(H_nogap_r, _MM_MAX_EPI8(pairscores_std,pairscores_alt));
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	debug15(print_vector_8(pairscores_std,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, _mm_max_epi8(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo])); /* Has 128 added already */
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(_mm_max_epu8(pairscores_std, pairscores_alt), all_128);
	debug15(print_vector_8(pairscores_best,rlo,c,"pairscores_std"));
	H_nogap_r = _mm_adds_epi8(H_nogap_r, pairscores_best);
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPGT_EPI8(E_r_gap,H_nogap_r); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_8(dir_horiz,rlo,c,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_r = _MM_MAX_EPI8(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_r = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_r, all_128), _mm_add_epi8(E_r_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif

	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi8(E_mask,15);
	E_mask = _mm256_slli_si256(E_mask,ONE_CHAR);
	E_mask = _mm256_insert_epi8(E_mask,save,16);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
#endif
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8_upper\n");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		   revp,uband,/*upperp*/true);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,uband,/*upperp*/true);
#endif
  
#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_8_upper_nonavx2(&directions_nogap_std,&directions_Egap_std,
					    this,rsequence,gsequence,gsequence_alt,
					    rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
					    open,extend,uband,jump_late_p,revp);

#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,/*lband*/0,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_8,
				/*upperp*/true,/*lowerp*/false);
#endif

#ifdef DEBUG2
  printf("Banded %s\n",revp ? "rev" : "fwd");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,revp,uband,/*upperp*/true);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,revp,uband,/*upperp*/true);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix8_compare_upper(matrix,matrix_std,rlength,glength,uband,
			       rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions8_compare_nogap_upper(matrix,*directions_nogap,directions_nogap_std,rlength,glength,uband);

  banded_directions8_compare_Egap_upper(matrix,*directions_Egap,directions_Egap_std,rlength,glength,uband);
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef DEBUG_AVX2
/* Designed for computation below the main diagonal, so no F loop or bottom masking needed */
/* Operates by rows */
Score8_T **
Dynprog_simd_8_lower_nonavx2 (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
			      T this, char *rsequence, char *gsequence, char *gsequence_alt,
			      int rlength, int glength,
			      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			      Mismatchtype_T mismatchtype, int open, int extend,
			      int lband, bool jump_late_p, bool revp) {
  Score8_T **matrix, *score_column;
  __m128i pairscores_std;
#ifdef HAVE_SSE4_1
  __m128i E_infinity;
#else
  __m128i pairscores_best, all_128, E_infinity_plus_128;
#endif
  __m128i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_vert;
  int glength_ceil, r, c;
  int clo, chigh;
  int na1, na2, na2_alt;
  Score8_T *pairscores[5], *pairscores_ptr;
  Pairdistance_T **pairdistance_array_type, score1, score2;


  debug2(printf("Dynprog_simd_8_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  glength_ceil = (int) ((glength + SIMD_NCHARS_NONAVX2)/SIMD_NCHARS_NONAVX2) * SIMD_NCHARS_NONAVX2;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("compute_scores_simd_8_byrows (lower): "));
  debug(printf("Lengths are %d and %d, so band is %d on left\n",rlength,glength,lband));
  debug(printf("Genome length rounded up to %d\n",glength_ceil));

  matrix = aligned_score8_alloc(glength_ceil,rlength,
				this->aligned_std.two.lower_matrix_ptrs,this->aligned_std.two.lower_matrix_space);
  *directions_nogap = aligned_directions8_alloc(glength_ceil,rlength,
						this->aligned_std.two.lower_directions_ptrs_0,this->aligned_std.two.lower_directions_space_0);
  *directions_Egap = aligned_directions8_alloc(glength_ceil,rlength,
					       this->aligned_std.two.lower_directions_ptrs_1,this->aligned_std.two.lower_directions_space_1);

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Egap)[r][0] = VERT;
    (*directions_nogap)[r][0] = VERT;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= SIMD_NCHARS_NONAVX2 && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_nogap)[0][c] = HORIZ;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),16);
  pairscores[1] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),16);
  pairscores[2] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),16);
  pairscores[3] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),16);
  pairscores[4] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,glength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
  c = 0; na2 = na2_alt = 4; /* 'N' */
#ifdef HAVE_SSE4_1
  pairscores[0][c] = (Score8_T) pairdistance_array_type[(int) 'A'][na2];
  pairscores[1][c] = (Score8_T) pairdistance_array_type[(int) 'C'][na2];
  pairscores[2][c] = (Score8_T) pairdistance_array_type[(int) 'G'][na2];
  pairscores[3][c] = (Score8_T) pairdistance_array_type[(int) 'T'][na2];
  pairscores[4][c] = (Score8_T) pairdistance_array_type[(int) 'N'][na2];
#else
  pairscores[0][c] = (Score8_T) pairdistance_array_type[(int) 'A'][na2] - 128;
  pairscores[1][c] = (Score8_T) pairdistance_array_type[(int) 'C'][na2] - 128;
  pairscores[2][c] = (Score8_T) pairdistance_array_type[(int) 'G'][na2] - 128;
  pairscores[3][c] = (Score8_T) pairdistance_array_type[(int) 'T'][na2] - 128;
  pairscores[4][c] = (Score8_T) pairdistance_array_type[(int) 'N'][na2] - 128;
#endif

  if (revp == false) {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[c-1];
      na2_alt = gsequence_alt[c-1];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score8_T) (score1 > score2) ? score1 : score2;
    }
  } else {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[1-c];
      na2_alt = gsequence_alt[1-c];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score8_T) (score1 > score2) ? score1 : score2;
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][c]),0,(glength_ceil-c)*sizeof(Score8_T));
#endif

  complement_dummy = _mm_set1_epi8(-1);

  gap_open = _mm_set1_epi8((Score8_T) open);
  gap_extend = _mm_set1_epi8((Score8_T) extend);

#ifdef HAVE_SSE4_1
  E_infinity = _mm_set1_epi8(POS_INFINITY_8);
#else
  E_infinity_plus_128 = _mm_set1_epi8(POS_INFINITY_8+128);
#endif
  if (jump_late_p) {
    for (clo = 0; clo <= glength; clo += SIMD_NCHARS_NONAVX2) {
      if ((chigh = clo + SIMD_NCHARS_NONAVX2 - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E >= H.  To fill in first row of each
	 column block with non-diags, make E == H. */
      E_mask = _mm_set1_epi8(1);

      E_c_gap = _mm_set1_epi8(NEG_INFINITY_8);
      H_nogap_c = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[rsequence[1-r]] : nt_to_int_array[rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the clo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	}

	debug15(print_vector_8(E_mask,clo,r,"E_mask"));
#ifdef HAVE_SSE4_1
	E_c_gap = _mm_min_epi8(E_c_gap,_mm_add_epi8(E_mask,E_infinity));
#else
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_c, gap_open);
	dir_vert = _mm_cmplt_epi8(E_c_gap,T1); /* E < H */
	dir_vert = _mm_andnot_si128(dir_vert,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
	debug15(print_vector_8(T1,clo,r,"T1"));
	debug15(print_vector_8(dir_vert,clo,r,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_c_gap = _mm_max_epi8(E_c_gap, T1); /* Compare H + open with horiz */
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _mm_min_epi8(E_c_gap,_mm_add_epi8(E_mask,E_infinity));
#elif 1
	E_c_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_c_gap = _mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128));
	E_c_gap = _mm_add_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(E_c_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_c,LAST_CHAR_NONAVX2);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_CHAR);
	H_nogap_c = _mm_or_si128(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(pairscores_std, all_128);
	debug15(print_vector_8(pairscores_best,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_best);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_c,clo,r,"H"));

	dir_vert = _mm_cmplt_epi8(E_c_gap,H_nogap_c); /* E < H */
	dir_vert = _mm_andnot_si128(dir_vert,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
	debug15(print_vector_8(dir_vert,clo,r,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_c = _mm_max_epi8(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_c = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_c, all_128), _mm_add_epi8(E_c_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c store"));
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
      }
    }

  } else {
    /* jump early */
    for (clo = 0; clo <= glength; clo += SIMD_NCHARS_NONAVX2) {
      if ((chigh = clo + SIMD_NCHARS_NONAVX2 - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E > H.  To fill in first row of each
	 column block with non-diags, make E > H. */
      E_mask = _mm_set1_epi8(1);

      E_c_gap = _mm_set1_epi8(NEG_INFINITY_8+1);
      H_nogap_c = _mm_set1_epi8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[rsequence[1-r]] : nt_to_int_array[rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _mm_set1_epi8(0);
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi8(0);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the clo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
#endif
	} else {
	  /* second or greater block of 8 */
	  X_prev_nogap = _mm_set1_epi8(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_NONAVX2);
	}

	debug15(print_vector_8(E_mask,clo,r,"E_mask"));
#ifdef HAVE_SSE4_1
	E_c_gap = _mm_min_epi8(E_c_gap,_mm_add_epi8(E_mask,E_infinity));
#else
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _mm_adds_epi8(H_nogap_c, gap_open);
	dir_vert = _mm_cmpgt_epi8(E_c_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
	debug15(print_vector_8(T1,clo,r,"T1"));
	debug15(print_vector_8(dir_vert,clo,r,"dir_Egap"));

	/* Compare H + open with vert */
#ifdef HAVE_SSE4_1
	E_c_gap = _mm_max_epi8(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_min_epi8(E_c_gap,_mm_add_epi8(E_mask,E_infinity));
#elif 1
	E_c_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_c_gap = _mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128));
	E_c_gap = _mm_add_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(E_c_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_c,LAST_CHAR_NONAVX2);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_CHAR);
	H_nogap_c = _mm_or_si128(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
#ifdef HAVE_SSE4_1
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(pairscores_std, all_128);
	debug15(print_vector_8(pairscores_best,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_best);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_c,clo,r,"H"));

	dir_vert = _mm_cmpgt_epi8(E_c_gap,H_nogap_c); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
	debug15(print_vector_8(dir_vert,clo,r,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_c = _mm_max_epi8(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_c = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_c, all_128), _mm_add_epi8(E_c_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c store"));
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);

	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
      }
    }
  }

#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = VERT;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8_lower\n");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		   revp,lband,/*upperp*/false);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,lband,/*upperp*/false);
#endif
  
#ifdef CHECK1
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Egap)[r][0] != DIAG);
    assert((*directions_nogap)[r][0] != DIAG);
  }
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef HAVE_SSE2
/* Designed for computation below the main diagonal, so no F loop or bottom masking needed */
/* Operates by rows */
Score8_T **
Dynprog_simd_8_lower (Direction8_T ***directions_nogap, Direction8_T ***directions_Egap,
		      T this, char *rsequence, char *gsequence, char *gsequence_alt,
		      int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		      int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		      Mismatchtype_T mismatchtype, int open, int extend,
		      int lband, bool jump_late_p, bool revp) {
  Score8_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i E_infinity;
#elif defined(HAVE_SSE4_1)
  __m128i E_infinity;
#else
  __m128i pairscores_best, all_128, E_infinity_plus_128;
#endif
#ifdef HAVE_AVX2
  __m256i pairscores_std;
  __m256i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_vert;
  Score8_T save;
#else
  __m128i pairscores_std;
  __m128i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_vert;
#endif
  int glength_ceil, r, c;
  int clo, chigh;
  int na1, na2, na2_alt;
  Score8_T *pairscores[5], *pairscores_ptr;
  Pairdistance_T **pairdistance_array_type, score1, score2;

#ifdef DEBUG_AVX2
  Score8_T **matrix_std;
  Direction8_T **directions_nogap_std, **directions_Egap_std;
  char na2_single;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif


  debug2(printf("Dynprog_simd_8_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_8_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  glength_ceil = (int) ((glength + SIMD_NCHARS)/SIMD_NCHARS) * SIMD_NCHARS;

#ifdef HAVE_SSE4_1
  pairdistance_array_type = pairdistance_array[mismatchtype];
#else
  /* Needed to use _mm_max_epu8 and _mm_min_epu8, instead of signed versions */
  pairdistance_array_type = pairdistance_array_plus_128[mismatchtype];
  all_128 = _mm_set1_epi8(128);
#endif
  
  debug(printf("compute_scores_simd_8_byrows (lower): "));
  debug(printf("Lengths are %d and %d, so band is %d on left\n",rlength,glength,lband));
  debug(printf("Genome length rounded up to %d\n",glength_ceil));

  matrix = aligned_score8_alloc(glength_ceil,rlength,
				this->aligned.two.lower_matrix_ptrs,this->aligned.two.lower_matrix_space);
  *directions_nogap = aligned_directions8_alloc(glength_ceil,rlength,
						this->aligned.two.lower_directions_ptrs_0,this->aligned.two.lower_directions_space_0);
  *directions_Egap = aligned_directions8_alloc(glength_ceil,rlength,
					       this->aligned.two.lower_directions_ptrs_1,this->aligned.two.lower_directions_space_1);

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Egap)[r][0] = VERT;
    (*directions_nogap)[r][0] = VERT;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= SIMD_NCHARS && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_nogap)[0][c] = HORIZ;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[1] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[2] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[3] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),ALIGN_SIZE);
  pairscores[4] = (Score8_T *) _mm_malloc(glength_ceil * sizeof(Score8_T),ALIGN_SIZE);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[1],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[2],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[3],0,glength_ceil*sizeof(Score8_T));
  memset((void *) pairscores[4],0,glength_ceil*sizeof(Score8_T));
#endif

  /* For non-SSE4.1, addition of 128 taken care of by using pairdistance_array_plus_128 above */
  c = 0; na2 = na2_alt = 4; /* 'N' */
#ifdef HAVE_SSE4_1
  pairscores[0][c] = (Score8_T) pairdistance_array_type[(int) 'A'][na2];
  pairscores[1][c] = (Score8_T) pairdistance_array_type[(int) 'C'][na2];
  pairscores[2][c] = (Score8_T) pairdistance_array_type[(int) 'G'][na2];
  pairscores[3][c] = (Score8_T) pairdistance_array_type[(int) 'T'][na2];
  pairscores[4][c] = (Score8_T) pairdistance_array_type[(int) 'N'][na2];
#else
  pairscores[0][c] = (Score8_T) pairdistance_array_type[(int) 'A'][na2] - 128;
  pairscores[1][c] = (Score8_T) pairdistance_array_type[(int) 'C'][na2] - 128;
  pairscores[2][c] = (Score8_T) pairdistance_array_type[(int) 'G'][na2] - 128;
  pairscores[3][c] = (Score8_T) pairdistance_array_type[(int) 'T'][na2] - 128;
  pairscores[4][c] = (Score8_T) pairdistance_array_type[(int) 'N'][na2] - 128;
#endif

  if (revp == false) {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[c-1];
      na2_alt = gsequence_alt[c-1];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score8_T) (score1 > score2) ? score1 : score2;
    }
  } else {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[1-c];
      na2_alt = gsequence_alt[1-c];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score8_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score8_T) (score1 > score2) ? score1 : score2;
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[1][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[2][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[3][c]),0,(glength_ceil-c)*sizeof(Score8_T));
  memset((void *) &(pairscores[4][c]),0,(glength_ceil-c)*sizeof(Score8_T));
#endif

  complement_dummy = _MM_SET1_EPI8(-1);

  gap_open = _MM_SET1_EPI8((Score8_T) open);
  gap_extend = _MM_SET1_EPI8((Score8_T) extend);

#ifdef HAVE_SSE4_1
  E_infinity = _MM_SET1_EPI8(POS_INFINITY_8);
#else
  E_infinity_plus_128 = _mm_set1_epi8(POS_INFINITY_8+128);
#endif
  if (jump_late_p) {
    for (clo = 0; clo <= glength; clo += SIMD_NCHARS) {
      if ((chigh = clo + SIMD_NCHARS - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E >= H.  To fill in first row of each
	 column block with non-diags, make E == H. */
      E_mask = _MM_SET1_EPI8(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_c_gap = _MM_SET1_EPI8(NEG_INFINITY_8);
      H_nogap_c = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[(int) rsequence[1-r]] : nt_to_int_array[(int) rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the clo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	} else {
	  /* second or greater block of 8 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[r-1][clo-1],LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	}

	debug15(print_vector_8(E_mask,clo,r,"E_mask"));
#ifdef HAVE_SSE4_1
	E_c_gap = _MM_MIN_EPI8(E_c_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#else
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_c, gap_open);
	dir_vert = _MM_CMPLT_EPI8(E_c_gap,T1); /* E < H */
	dir_vert = _MM_ANDNOT_SI(dir_vert,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_8(T1,clo,r,"T1"));
	debug15(print_vector_8(dir_vert,clo,r,"dir_Egap"));

#ifdef HAVE_SSE4_1
	E_c_gap = _MM_MAX_EPI8(E_c_gap, T1); /* Compare H + open with horiz */
	E_c_gap = _MM_ADDS_EPI8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _MM_MIN_EPI8(E_c_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#elif 1
	E_c_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_c_gap = _mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128));
	E_c_gap = _mm_add_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (horiz + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(E_c_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_c,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_c,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_c = _mm256_slli_si256(H_nogap_c,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_c,LAST_CHAR_SHIFT);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_CHAR);
#endif
	H_nogap_c = _MM_OR_SI(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _MM_ADDS_EPI8(H_nogap_c, pairscores_std);
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(pairscores_std, all_128);
	debug15(print_vector_8(pairscores_best,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_best);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_c,clo,r,"H"));

	dir_vert = _MM_CMPLT_EPI8(E_c_gap,H_nogap_c); /* E < H */
	dir_vert = _MM_ANDNOT_SI(dir_vert,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_8(dir_vert,clo,r,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_c = _MM_MAX_EPI8(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_c = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_c, all_128), _mm_add_epi8(E_c_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[clo]), H_nogap_c);
#else
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);
#endif


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi8(E_mask,15);
	E_mask = _mm256_slli_si256(E_mask,ONE_CHAR);
	E_mask = _mm256_insert_epi8(E_mask,save,16);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
#endif
      }
    }

  } else {
    /* jump early */
    for (clo = 0; clo <= glength; clo += SIMD_NCHARS) {
      if ((chigh = clo + SIMD_NCHARS - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E > H.  To fill in first row of each
	 column block with non-diags, make E > H. */
      E_mask = _MM_SET1_EPI8(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_c_gap = _MM_SET1_EPI8(NEG_INFINITY_8+1);
      H_nogap_c = _MM_SET1_EPI8(NEG_INFINITY_8-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[(int) rsequence[1-r]] : nt_to_int_array[(int) rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),NEG_INFINITY_8,LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(NEG_INFINITY_8); /* works if we start outside the clo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	} else {
	  /* second or greater block of 8 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi8(_mm256_setzero_si256(),matrix[r-1][clo-1],LAST_CHAR_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi8(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_CHAR_SHIFT);
#endif
	}

	debug15(print_vector_8(E_mask,clo,r,"E_mask"));
#ifdef HAVE_SSE4_1
	E_c_gap = _MM_MIN_EPI8(E_c_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#else
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI8(H_nogap_c, gap_open);
	dir_vert = _MM_CMPGT_EPI8(E_c_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_8(T1,clo,r,"T1"));
	debug15(print_vector_8(dir_vert,clo,r,"dir_Egap"));

	/* Compare H + open with vert */
#ifdef HAVE_SSE4_1
	E_c_gap = _MM_MAX_EPI8(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _MM_ADDS_EPI8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _MM_MIN_EPI8(E_c_gap,_MM_ADD_EPI8(E_mask,E_infinity));
#elif 1
	E_c_gap = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128)), all_128);
	E_c_gap = _mm_adds_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(_mm_add_epi8(E_c_gap, all_128),_mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#else
	/* Try to avoid unnecessary shifts by 128, but overflows */
	E_c_gap = _mm_max_epu8(_mm_add_epi8(E_c_gap, all_128), _mm_add_epi8(T1, all_128));
	E_c_gap = _mm_add_epi8(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_sub_epi8(_mm_min_epu8(E_c_gap, _mm_add_epi8(E_mask,E_infinity_plus_128)), all_128);
#endif
	debug15(print_vector_8(E_c_gap,clo,r,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi8(_mm256_setzero_si256(),_mm256_extract_epi8(H_nogap_c,SIMD_NCHARS-1),LAST_CHAR_INSERT);
	X_prev_nogap = _mm256_insert_epi8(X_prev_nogap,_mm256_extract_epi8(H_nogap_c,MID_CHAR_INSERT-1),MID_CHAR_INSERT);
	H_nogap_c = _mm256_slli_si256(H_nogap_c,ONE_CHAR);
#else
	T1 = _mm_srli_si128(H_nogap_c,LAST_CHAR_SHIFT);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_CHAR);
#endif
	H_nogap_c = _MM_OR_SI(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _MM_ADDS_EPI8(H_nogap_c, pairscores_std);
#elif defined(HAVE_SSE4_1)
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	debug15(print_vector_8(pairscores_std,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo])); /* Has 128 added already */
	pairscores_best = _mm_sub_epi8(pairscores_std, all_128);
	debug15(print_vector_8(pairscores_best,clo,r,"pairscores_std"));
	H_nogap_c = _mm_adds_epi8(H_nogap_c, pairscores_best);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_8(H_nogap_c,clo,r,"H"));

	dir_vert = _MM_CMPGT_EPI8(E_c_gap,H_nogap_c); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_8(dir_vert,clo,r,"dir_nogap"));


#ifdef HAVE_SSE4_1
	H_nogap_c = _MM_MAX_EPI8(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
#else
	/* Compare H + pairscores with horiz + extend */
	H_nogap_c = _mm_sub_epi8(_mm_max_epu8(_mm_add_epi8(H_nogap_c, all_128), _mm_add_epi8(E_c_gap, all_128)), all_128);
#endif
	debug15(print_vector_8(H_nogap_c,clo,r,"H_nogap_c store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[clo]), H_nogap_c);
#else
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);
#endif

	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi8(E_mask,15);
	E_mask = _mm256_slli_si256(E_mask,ONE_CHAR);
	E_mask = _mm256_insert_epi8(E_mask,save,16);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_CHAR);
#endif
      }
    }
  }

#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = VERT;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_8_lower\n");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		   revp,lband,/*upperp*/false);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,lband,/*upperp*/false);
#endif
  
#ifdef CHECK1
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Egap)[r][0] != DIAG);
    assert((*directions_nogap)[r][0] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_8_lower_nonavx2(&directions_nogap_std,&directions_Egap_std,
					    this,rsequence,gsequence,gsequence_alt,
					    rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
					    open,extend,lband,jump_late_p,revp);
#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,lband,/*uband*/0,jump_late_p,revp,/*saturation*/NEG_INFINITY_8,
				/*upperp*/false,/*lowerp*/true);
#endif

#ifdef DEBUG2
  printf("Banded %s\n",revp ? "rev" : "fwd");
  Matrix8_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		   revp,lband,/*upperp*/false);
  Directions8_print_ud(*directions_nogap,*directions_Egap,
		       rlength,glength,rsequence,gsequence,gsequence_alt,
		       revp,lband,/*upperp*/false);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix8_compare_lower(matrix,matrix_std,rlength,glength,lband,
			       rsequence,gsequence,gsequence_alt,
			       goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions8_compare_nogap_lower(matrix,*directions_nogap,directions_nogap_std,rlength,glength,lband);

  banded_directions8_compare_Egap_lower(matrix,*directions_Egap,directions_Egap_std,rlength,glength,lband);
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef DEBUG_AVX2
Score16_T **
Dynprog_simd_16_nonavx2 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
			 Direction16_T ***directions_Fgap,
			 T this, char *rsequence, char *gsequence, char *gsequence_alt,
			 int rlength, int glength,
			 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			 Mismatchtype_T mismatchtype, int open, int extend,
			 int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF; /* Need to have the ability to go past NEG_INFINITY */
  Score16_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, extend_ladder, extend_chunk, complement_dummy;
  __m128i dir_horiz;
  __m128i ramp, ramp_chunk, lband_vector, filter, ones;
  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score16_T *pairscores_col0;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr, pairscore, pairscore0;
  Pairdistance_T **pairdistance_array_type;


  debug2(printf("Dynprog_simd_16.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16.  jump_late_p %d, open %d, extend %d, lband %d, uband %d\n",
		 jump_late_p,open,extend,lband,uband));

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS_NONAVX2)/SIMD_NSHORTS_NONAVX2) * SIMD_NSHORTS_NONAVX2;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned_std.one.matrix_ptrs,this->aligned_std.one.matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned_std.one.directions_ptrs_0,this->aligned_std.one.directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned_std.one.directions_ptrs_1,this->aligned_std.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions16_calloc(rlength_ceil,glength,
						 this->aligned_std.one.directions_ptrs_2,this->aligned_std.one.directions_space_2);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NSHORTS_NONAVX2 && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);

  lband_ceil = (int) ((lband + SIMD_NSHORTS_NONAVX2)/SIMD_NSHORTS_NONAVX2) * SIMD_NSHORTS_NONAVX2;
  pairscores_col0 = (Score16_T *) _mm_malloc(lband_ceil * sizeof(Score16_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score16_T));
#endif


  pairscores_col0[0] = (Score16_T) 0;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score16_T) 0;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score16_T) NEG_INFINITY_16;
  }
#endif

  r = 0; na1 = 'N';
  pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
#endif

  complement_dummy = _mm_set1_epi16(-1);
  
  FF = (int *) MALLOCA((glength + 1) * sizeof(int));

  gap_open = _mm_set1_epi16((Score16_T) open);
  gap_extend = _mm_set1_epi16((Score16_T) extend);


#ifndef INFINITE_INITIAL_GAP_PENALTY
  lband_vector = _mm_set1_epi16(lband);
  ramp = _mm_setr_epi16(1,2,3,4,5,6,7,8);
  extend_ladder = _mm_setr_epi16(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend);
  ramp_chunk = _mm_set1_epi16(SIMD_NSHORTS_NONAVX2);
  extend_chunk = _mm_set1_epi16(SIMD_NSHORTS_NONAVX2*extend);
#endif

  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NSHORTS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16-open);
	filter = _mm_cmpgt_epi16(ramp,lband_vector);
	H_nogap_r = _mm_and_si128(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _mm_adds_epi16(E_r_gap,gap_open);

	ramp = _mm_adds_epi16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi16(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E >= H.  To fill in first column of each
	   row block with non-diags, make E == H. */
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16);
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   equal to E. */
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_8-open); /* To compensate for T1 = H + open */
	filter = _mm_cmpgt_epi16(ramp,lband_vector);
	H_nogap_r = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),_mm_and_si128(filter,H_nogap_r));
	E_r_gap = _mm_adds_epi16(H_nogap_r,gap_open);
	ramp = _mm_adds_epi16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi16(extend_ladder,extend_chunk); /* Prepare for next block */
#endif

      } else {
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16);
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_set1_epi16(0);
	  } else {
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _mm_set1_epi16(0);
#else
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	  }
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi16(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) <= c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
	    /* No need to fix for non-SSE4.1: pairscore -= 128; */
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_16) {
	      score_column[rhigh_calc] = NEG_INFINITY_16; /* Saturation */
	    } else if (score > POS_INFINITY_16) {
	      score_column[rhigh_calc] = POS_INFINITY_16; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score16_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }

  } else {
    /* jump early */
#if defined(ZERO_INITIAL_GAP_PENALTY) || defined(INFINITE_INITIAL_GAP_PENALTY)
    /* No need for ones */
#else
    ones = _mm_set1_epi16(1);
#endif
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NSHORTS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16-open);
	filter = _mm_cmpgt_epi16(ramp,lband_vector);
	H_nogap_r = _mm_and_si128(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _mm_adds_epi16(E_r_gap,gap_open);

	ramp = _mm_adds_epi16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi16(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E > H.  To fill in first column of each
	   row block with non-diags, make E > H. */
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16+1);	
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   less than E. */
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open+1); /* To compensate for T1 = H + open */
	filter = _mm_cmpgt_epi16(ramp,lband_vector);
	H_nogap_r = _mm_or_si128(_mm_andnot_si128(filter,extend_ladder),_mm_and_si128(filter,H_nogap_r));
	E_r_gap = _mm_adds_epi16(H_nogap_r,gap_open);
	H_nogap_r = _mm_subs_epi16(H_nogap_r,ones);    /* To ensure H < E */
	ramp = _mm_adds_epi16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _mm_adds_epi16(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
	
      } else {
	E_r_gap = _mm_set1_epi16(NEG_INFINITY_16+1);
	H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_set1_epi16(0);
	  } else {
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _mm_set1_epi16(0);
#else
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	  }
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi16(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi16(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
	    /* No need to fix for non-SSE4.1: pairscore -= 128; */
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_16) {
	      score_column[rhigh_calc] = NEG_INFINITY_16; /* Saturation */
	    } else if (score > POS_INFINITY_16) {
	      score_column[rhigh_calc] = POS_INFINITY_16; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score16_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
  (*directions_Fgap)[0][1] = VERT;
#endif


#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16\n");
  Matrix16_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		 revp,lband,uband);
  Directions16_print(*directions_nogap,*directions_Egap,*directions_Fgap,
		     rlength,glength,rsequence,gsequence,gsequence_alt,revp,lband,uband);
#endif

#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Fgap)[0][r] != DIAG);
    assert((*directions_nogap)[0][r] != DIAG);
  }
#endif

  FREEA(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
  }
#endif




#if defined(HAVE_SSE2)
/* Modified from Dynprog_simd_16_upper.  Operates by columns. */
Score16_T **
Dynprog_simd_16 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
		 Direction16_T ***directions_Fgap,
		 T this, char *rsequence, char *gsequence, char *gsequence_alt,
		 int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		 int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		 Mismatchtype_T mismatchtype, int open, int extend,
		 int lband, int uband, bool jump_late_p, bool revp) {
  int c_gap, last_nogap, score, *FF; /* Need to have the ability to go past NEG_INFINITY */
  Score16_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i pairscores_std, pairscores_alt;
  __m256i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_horiz;
#if defined(ZERO_INITIAL_GAP_PENALTY)
  __m256i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter;
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
#else
  __m256i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter, ones;
#endif

#else
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
#if defined(ZERO_INITIAL_GAP_PENALTY)
  __m128i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter;
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
#else
  __m128i ramp, ramp_chunk, extend_ladder, extend_chunk, lband_vector, filter, ones;
#endif

#endif

  int rlength_ceil, lband_ceil, r, c;
  int rlo, rhigh, rlo_calc, rhigh_calc;
  int na1, na2, na2_alt;
  Score16_T *pairscores_col0;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr, pairscore, pairscore0;
  Pairdistance_T **pairdistance_array_type;

#if defined(DEBUG_AVX2)
  Score16_T **matrix_std;
  Direction16_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
#endif


  debug2(printf("Dynprog_simd_16.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16.  jump_late_p %d, open %d, extend %d, lband %d, uband %d\n",
		 jump_late_p,open,extend,lband,uband));

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned.one.matrix_ptrs,this->aligned.one.matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned.one.directions_ptrs_0,this->aligned.one.directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned.one.directions_ptrs_1,this->aligned.one.directions_space_1);
  /* Need to calloc to save time in F loop */
  *directions_Fgap = aligned_directions16_calloc(rlength_ceil,glength,
						 this->aligned.one.directions_ptrs_2,this->aligned.one.directions_space_2);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NSHORTS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);

  lband_ceil = (int) ((lband + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairscores_col0 = (Score16_T *) _mm_malloc(lband_ceil * sizeof(Score16_T),ALIGN_SIZE);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score16_T));
#endif


  pairscores_col0[0] = (Score16_T) 0;
  /* Initialization just to lband causes errors in dir_horiz for Egap */
#ifdef ZERO_INITIAL_GAP_PENALTY
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score16_T) 0;
  }
#else
  for (r = 1; r < lband_ceil; r++) {
    pairscores_col0[r] = (Score16_T) NEG_INFINITY_16;
  }
#endif

  r = 0; na1 = 'N';
  pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
#endif

  complement_dummy = _MM_SET1_EPI16(-1);
  
  FF = (int *) MALLOCA((glength + 1) * sizeof(int));

  gap_open = _MM_SET1_EPI16((Score16_T) open);
  gap_extend = _MM_SET1_EPI16((Score16_T) extend);

#ifndef INFINITE_INITIAL_GAP_PENALTY
#ifdef HAVE_AVX2
  ramp = _mm256_setr_epi16(1,2,3,4,5,6,7,8);
  extend_ladder = _mm256_setr_epi16(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend);
#else
  ramp = _mm_setr_epi16(1,2,3,4,5,6,7,8);
  extend_ladder = _mm_setr_epi16(extend,2*extend,3*extend,4*extend,5*extend,6*extend,7*extend,8*extend);
#endif
  lband_vector = _MM_SET1_EPI16(lband);
  ramp_chunk = _MM_SET1_EPI16(SIMD_NSHORTS);
  extend_chunk = _MM_SET1_EPI16(SIMD_NSHORTS*extend);
#endif

  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16-open);
	filter = _MM_CMPGT_EPI16(ramp,lband_vector);
	H_nogap_r = _MM_AND_SI(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _MM_OR_SI(_mm_andnot_si128(filter,extend_ladder),H_nogap_r);
	E_r_gap = _MM_ADDS_EPI16(E_r_gap,gap_open);

	ramp = _MM_ADDS_EPI16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI16(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E >= H.  To fill in first column of each
	   row block with non-diags, make E == H. */
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16);
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   equal to E. */
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_8-open); /* To compensate for T1 = H + open */
	filter = _MM_CMPGT_EPI16(ramp,lband_vector);
	H_nogap_r = _MM_OR_SI(_mm_andnot_si128(filter,extend_ladder),_mm_and_si128(filter,H_nogap_r));
	E_r_gap = _MM_ADDS_EPI16(H_nogap_r,gap_open);
	ramp = _MM_ADDS_EPI16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI16(extend_ladder,extend_chunk); /* Prepare for next block */
#endif

      } else {
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16);
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  if (rlo == 0) {
	    X_prev_nogap = _mm256_setzero_si256();
	  } else {
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
	  }
#else
	  if (rlo == 0) {
	    X_prev_nogap = _mm_setzero_si128();
	  } else {
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

#ifdef HAVE_AVX2
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_SHORT_INSERT);
	  }

#else
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
	  }
#endif
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPLT_EPI16(E_r_gap,T1); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _MM_MAX_EPI16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_r,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_r,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi16(H_nogap_r, _mm256_max_epi16(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPLT_EPI16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _MM_MAX_EPI16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif


	/* F loop */
	if ((rlo_calc = rlo) <= c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
	    /* No need to fix for non-SSE4.1: pairscore -= 128; */
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_16) {
	      score_column[rhigh_calc] = NEG_INFINITY_16; /* Saturation */
	    } else if (score > POS_INFINITY_16) {
	      score_column[rhigh_calc] = POS_INFINITY_16; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score16_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  debug3(printf("At c %d, uband %d, reading c_gap %d\n",c,uband,FF[c]));
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */

	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ >= (score = last_nogap + open /* + extend */)) {  /* Use >= for jump late */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap >= last_nogap) {  /* Use >= for jump late */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
#ifdef HAVE_AVX2
	H_nogap_r = _mm256_load_si256((__m256i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#else
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#endif
      }
    }

  } else {
    /* jump early */
#if defined(ZERO_INITIAL_GAP_PENALTY) || defined(INFINITE_INITIAL_GAP_PENALTY)
    /* No need for ones */
#else
    ones = _MM_SET1_EPI16(1);
#endif
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      if ((c = rlo - lband) < 0) {
	c = 0;

#if defined(ZERO_INITIAL_GAP_PENALTY)
	/* Initial H in column 0 determined by zeroed out H.  E needs to equal gap_open for column 1. */
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16-open);
	filter = _MM_CMPGT_EPI16(ramp,lband_vector);
	H_nogap_r = _MM_AND_SI(filter,E_r_gap); /* Use zeros for score */

	E_r_gap = _MM_OR_SI(_MM_ANDNOT_SI(filter,extend_ladder),H_nogap_r);
	E_r_gap = _MM_ADDS_EPI16(E_r_gap,gap_open);

	ramp = _MM_ADDS_EPI16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI16(extend_ladder,extend_chunk); /* Prepare for next block */
#elif defined(INFINITE_INITIAL_GAP_PENALTY)
	/* dir_horiz tests if E > H.  To fill in first column of each
	   row block with non-diags, make E > H. */
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16+1);	
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
#else
	/* Initial H in column 0 determined by E, which needs to equal
	   gap_open + extend_ladder for column 1.  H is free to be set
	   less than E. */
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open+1); /* To compensate for T1 = H + open */
	filter = _MM_CMPGT_EPI16(ramp,lband_vector);
	H_nogap_r = _MM_OR_SI(_MM_ANDNOT_SI(filter,extend_ladder),_MM_AND_SI(filter,H_nogap_r));
	E_r_gap = _MM_ADDS_EPI16(H_nogap_r,gap_open);
	H_nogap_r = _MM_SUBS_EPI16(H_nogap_r,ones);    /* To ensure H < E */
	ramp = _MM_ADDS_EPI16(ramp,ramp_chunk); /* Prepare for next block */
	extend_ladder = _MM_ADDS_EPI16(extend_ladder,extend_chunk); /* Prepare for next block */
#endif
	
      } else {
	E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16+1);
	H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */
      }

      for ( ; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  pairscores_std_ptr = pairscores_alt_ptr = pairscores_col0;

#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  if (rlo == 0) {
            X_prev_nogap = _mm256_setzero_si256();
	  } else {
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
	  }
#else
	  if (rlo == 0) {
            X_prev_nogap = _mm_setzero_si128();
	  } else {
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
	  }
#endif

	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	  pairscores_std_ptr = pairscores[na2];
	  pairscores_alt_ptr = pairscores[na2_alt];

#ifdef HAVE_AVX2
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_SHORT_INSERT);
	  }

#else
	  if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	    X_prev_nogap = _MM_SETZERO_SI();
#else
	    X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16);
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	  } else {
	    /* second or greater block of 16 */
	    X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	    X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
	  }
#endif
	}

	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPGT_EPI16(E_r_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _MM_MAX_EPI16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_r,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_r,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi16(H_nogap_r, _mm256_max_epi16(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPGT_EPI16(E_r_gap,H_nogap_r); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _MM_MAX_EPI16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif


	/* F loop */
	if ((rlo_calc = rlo) < c - uband) {
	  rlo_calc = c - uband;
	}
	if ((rhigh_calc = rhigh) >= c + lband) {
	  rhigh_calc = c + lband;
	  if (c > 0) {
	    /* Set bottom values to DIAG (not HORIZ) to prevent going outside of lband */
	    pairscore = pairscores[na2][rhigh_calc];
	    if ((pairscore0 = pairscores[(int) na2_alt][rhigh_calc]) > pairscore) {
	      pairscore = pairscore0;
	    }
	    /* No need to fix for non-SSE4.1: pairscore -= 128; */
	    if ((score = (int) matrix[c-1][rhigh_calc-1] + (int) pairscore) < NEG_INFINITY_16) {
	      score_column[rhigh_calc] = NEG_INFINITY_16; /* Saturation */
	    } else if (score > POS_INFINITY_16) {
	      score_column[rhigh_calc] = POS_INFINITY_16; /* Saturation */
	    } else {
	      score_column[rhigh_calc] = (Score16_T) score;
	    }
	    (*directions_Egap)[c][rhigh_calc] = DIAG;
	    (*directions_nogap)[c][rhigh_calc] = DIAG;
	  }
	}

	debug3(printf("F loop: rlo %d, rhigh %d, c %d, lband %d, uband %d => rlo_calc %d, rhigh_calc %d\n",
		      rlo,rhigh,rlo_calc,c,lband,uband,rhigh_calc));

	if (rlo == 0) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else if (c >= rlo + uband) {
	  c_gap = NEG_INFINITY_INT;
	  last_nogap = NEG_INFINITY_INT;
	} else {
	  c_gap = FF[c];
	  last_nogap = (int) score_column[rlo_calc-1];
	}

	if ((r = rlo_calc) == c - uband) {
	  /* Handle top value as a special case to prevent going outside of uband */
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  score = last_nogap + open /* + extend */;
	  c_gap = score + extend;
	  /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  r++;
	}

	/* score_ptr = &(score_column[rlo_calc]); -- Also possible, but less transparent */
	for ( ; r <= rhigh_calc; r++) {
	  /* FGAP */
	  debug3(printf("Fgap at r %d, c %d: c_gap + extend %d vs last_nogap + open + extend %d\n",
			r,c,c_gap + extend,last_nogap + open + extend));
	  if (c_gap /* + extend */ > (score = last_nogap + open /* + extend */)) {  /* Use > for jump early */
	    c_gap += extend;
	    (*directions_Fgap)[c][r] = VERT;
	  } else {
	    c_gap = score + extend;
	    /* (*directions_Fgap)[c][r] = DIAG: -- Already initialized to DIAG */
	  }
	  
	  /* NOGAP */
	  last_nogap = (int) score_column[r];
	  debug3(printf("assign nogap at r %d, c %d: H/E %d vs vert + extend %d\n",r,c,last_nogap,c_gap));
	  if (c_gap > last_nogap) {  /* Use > for jump early */
	    last_nogap = c_gap;
	    score_column[r] = (c_gap < NEG_INFINITY_16) ? NEG_INFINITY_16 : (Score16_T) c_gap; /* Saturation */
	    (*directions_nogap)[c][r] = VERT;
	  }
	}

	FF[c] = c_gap;
	debug3(printf("At c %d, storing c_gap %d\n",c,FF[c]));
#ifdef HAVE_AVX2
	H_nogap_r = _mm256_load_si256((__m256i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#else
	H_nogap_r = _mm_load_si128((__m128i *) &(score_column[rlo])); /* Need to reload because of changes by F loop */
#endif
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
  (*directions_Fgap)[0][1] = VERT;
#endif


#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16\n");
  Matrix16_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		 revp,lband,uband);
  Directions16_print(*directions_nogap,*directions_Egap,*directions_Fgap,
			     rlength,glength,rsequence,gsequence,gsequence_alt,
			     revp,lband,uband);
#endif

#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Fgap)[0][r] != DIAG);
    assert((*directions_nogap)[0][r] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_16_nonavx2(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				       this,rsequence,gsequence,gsequence_alt,
				       rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				       open,extend,lband,uband,jump_late_p,revp);
#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,lband,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_16,
				/*upperp*/true,/*lowerp*/true);
#endif

#ifdef DEBUG2
  printf("Banded\n");
  Matrix16_print(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		 revp,lband,uband);
  Directions16_print(*directions_nogap,*directions_Egap,*directions_Fgap,
		     rlength,glength,rsequence,gsequence,gsequence_alt,revp,lband,uband);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix16_compare(matrix,matrix_std,rlength,glength,lband,uband,
			  rsequence,gsequence,gsequence_alt,
			  goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions16_compare_nogap(*directions_nogap,directions_nogap_std,rlength,glength,lband,uband);
  banded_directions16_compare_Egap(*directions_Egap,directions_Egap_std,rlength,glength,lband,uband);
  banded_directions16_compare_Fgap(*directions_Fgap,directions_Fgap_std,rlength,glength,lband,uband);
#endif

  FREEA(FF);
  _mm_free(pairscores_col0);
  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
  }
#endif


#ifdef DEBUG_AVX2
/* Designed for computation above the diagonal, so no F loop or bottom masking needed */
/* Operates by columns */
Score16_T **
Dynprog_simd_16_upper_nonavx2 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
			       T this, char *rsequence, char *gsequence, char *gsequence_alt,
			       int rlength, int glength,
			       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			       Mismatchtype_T mismatchtype, int open, int extend,
			       int uband, bool jump_late_p, bool revp) {
  Score16_T **matrix, *score_column;
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, E_infinity, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
  int rlength_ceil, r, c;
  int rlo, rhigh;
  int na1, na2, na2_alt;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;


  debug2(printf("Dynprog_simd_16_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS_NONAVX2)/SIMD_NSHORTS_NONAVX2) * SIMD_NSHORTS_NONAVX2;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned_std.two.upper_matrix_ptrs,this->aligned_std.two.upper_matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned_std.two.upper_directions_ptrs_0,this->aligned_std.two.upper_directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned_std.two.upper_directions_ptrs_1,this->aligned_std.two.upper_directions_space_1);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NSHORTS_NONAVX2 && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score16_T));
#endif

  r = 0; na1 = 'N';
  pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
#endif

  complement_dummy = _mm_set1_epi16(-1);
  
  gap_open = _mm_set1_epi16((Score16_T) open);
  gap_extend = _mm_set1_epi16((Score16_T) extend);

  E_infinity = _mm_set1_epi16(POS_INFINITY_16);
  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NSHORTS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H.  To fill in first column of each
	 row block with non-diags, could make E == H.  But irrelevant,
	 because these are above the diagonal. */
      E_mask = _mm_set1_epi16(1);

      E_r_gap = _mm_set1_epi16(NEG_INFINITY_16);
      H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	}

	debug15(print_vector_16(E_mask,rlo,c,"E_mask"));
	E_r_gap = _mm_min_epi16(E_r_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmplt_epi16(E_r_gap,T1); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_min_epi16(E_r_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmplt_epi16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _mm_andnot_si128(dir_horiz,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
      }
    }

  } else {
    /* jump early */
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS_NONAVX2) {
      if ((rhigh = rlo + SIMD_NSHORTS_NONAVX2 - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with non-diags, could make E > H.  But irrelevant,
	 because these are above the diagonal. */
      E_mask = _mm_set1_epi16(1);

      E_r_gap = _mm_set1_epi16(NEG_INFINITY_16+1);
      H_nogap_r = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[gsequence[1-c]] : nt_to_int_array[gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[gsequence_alt[1-c]] : nt_to_int_array[gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	}

	debug15(print_vector_16(E_mask,rlo,c,"E_mask"));
	E_r_gap = _mm_min_epi16(E_r_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_r, gap_open);
	dir_horiz = _mm_cmpgt_epi16(E_r_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _mm_max_epi16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _mm_adds_epi16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _mm_min_epi16(E_r_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_NONAVX2);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
	H_nogap_r = _mm_or_si128(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _mm_cmpgt_epi16(E_r_gap,H_nogap_r); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _mm_max_epi16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);


	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
      }
    }
  }

#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16_upper\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		    revp,uband,/*upperp*/true);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
			rlength,glength,rsequence,gsequence,gsequence_alt,
			revp,uband,/*upperp*/true);
#endif

#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef HAVE_SSE2
/* Designed for computation above the diagonal, so no F loop or bottom masking needed */
/* Operates by columns */
Score16_T **
Dynprog_simd_16_upper (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
		       T this, char *rsequence, char *gsequence, char *gsequence_alt,
		       int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		       Mismatchtype_T mismatchtype, int open, int extend,
		       int uband, bool jump_late_p, bool revp) {
  Score16_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i pairscores_std, pairscores_alt;
  __m256i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, E_infinity, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_horiz;
  Score16_T save;
#else
  __m128i pairscores_std, pairscores_alt;
  __m128i H_nogap_r, X_prev_nogap, E_r_gap, E_mask, E_infinity, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_horiz;
#endif
  int rlength_ceil, r, c;
  int rlo, rhigh;
  int na1, na2, na2_alt;
  Score16_T *pairscores[5], *pairscores_std_ptr, *pairscores_alt_ptr;
  Pairdistance_T **pairdistance_array_type;

#ifdef DEBUG_AVX2
  Score16_T **matrix_std;
  Direction16_T **directions_nogap_std, **directions_Egap_std;
  char na2_single;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif


  debug2(printf("Dynprog_simd_16_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16_upper.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  rlength_ceil = (int) ((rlength + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_bycols (upper): "));
  debug(printf("Lengths are %d and %d, so band is %d on right\n",rlength,glength,uband));
  debug(printf("Query length rounded up to %d\n",rlength_ceil));

  matrix = aligned_score16_alloc(rlength_ceil,glength,
				 this->aligned.two.upper_matrix_ptrs,this->aligned.two.upper_matrix_space);
  *directions_nogap = aligned_directions16_alloc(rlength_ceil,glength,
						 this->aligned.two.upper_directions_ptrs_0,this->aligned.two.upper_directions_space_0);
  *directions_Egap = aligned_directions16_alloc(rlength_ceil,glength,
						this->aligned.two.upper_directions_ptrs_1,this->aligned.two.upper_directions_space_1);

#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= uband && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_Egap)[c][0] = HORIZ;
    (*directions_nogap)[c][0] = HORIZ;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= SIMD_NSHORTS && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_nogap)[0][r] = VERT;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[1] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[2] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[3] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[4] = (Score16_T *) _mm_malloc(rlength_ceil * sizeof(Score16_T),ALIGN_SIZE);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,rlength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,rlength_ceil*sizeof(Score16_T));
#endif

  r = 0; na1 = 'N';
  pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
  pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
  pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
  pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
  pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];

  if (revp == false) {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[r-1];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  } else {
    for (r = 1; r <= rlength; r++) {
      na1 = (int) rsequence[1-r];
      pairscores[0][r] = (Score16_T) pairdistance_array_type[na1][(int) 'A'];
      pairscores[1][r] = (Score16_T) pairdistance_array_type[na1][(int) 'C'];
      pairscores[2][r] = (Score16_T) pairdistance_array_type[na1][(int) 'G'];
      pairscores[3][r] = (Score16_T) pairdistance_array_type[na1][(int) 'T'];
      pairscores[4][r] = (Score16_T) pairdistance_array_type[na1][(int) 'N'];
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][r]),0,(rlength_ceil-r)*sizeof(Score16_T));
#endif

  complement_dummy = _MM_SET1_EPI16(-1);
  
  gap_open = _MM_SET1_EPI16((Score16_T) open);
  gap_extend = _MM_SET1_EPI16((Score16_T) extend);

  E_infinity = _MM_SET1_EPI16(POS_INFINITY_16);
  if (jump_late_p) {
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E >= H.  To fill in first column of each
	 row block with non-diags, could make E == H.  But irrelevant,
	 because these are above the diagonal. */
      E_mask = _MM_SET1_EPI16(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16);
      H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	} else {
	  /* second or greater block of 16 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	}

	debug15(print_vector_16(E_mask,rlo,c,"E_mask"));
	E_r_gap = _MM_MIN_EPI16(E_r_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPLT_EPI16(E_r_gap,T1); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _MM_MAX_EPI16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _MM_MIN_EPI16(E_r_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_r,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_r,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi16(H_nogap_r, _mm256_max_epi16(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPLT_EPI16(E_r_gap,H_nogap_r); /* E < H */
	dir_horiz = _MM_ANDNOT_SI(dir_horiz,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _MM_MAX_EPI16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif

	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi16(E_mask,7);
	E_mask = _mm256_slli_si256(E_mask,ONE_SHORT);
	E_mask = _mm256_insert_epi16(E_mask,save,8);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
#endif
      }
    }

  } else {
    /* jump early */
    for (rlo = 0; rlo <= rlength; rlo += SIMD_NSHORTS) {
      if ((rhigh = rlo + SIMD_NSHORTS - 1) > rlength) {
	rhigh = rlength;
      }

      /* dir_horiz tests if E > H.  To fill in first column of each
	 row block with non-diags, could make E > H.  But irrelevant,
	 because these are above the diagonal. */
      E_mask = _MM_SET1_EPI16(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_r_gap = _MM_SET1_EPI16(NEG_INFINITY_16+1);
      H_nogap_r = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (c = rlo; c <= rhigh + uband && c <= glength; c++) {
	score_column = matrix[c];

	if (c == 0) {
	  na2 = na2_alt = 4; /* 'N' */
	} else {
	  na2 = revp ? nt_to_int_array[(int) gsequence[1-c]] : nt_to_int_array[(int) gsequence[c-1]];
	  na2_alt = revp ? nt_to_int_array[(int) gsequence_alt[1-c]] : nt_to_int_array[(int) gsequence_alt[c-1]];
	}
	pairscores_std_ptr = pairscores[na2];
	pairscores_alt_ptr = pairscores[na2_alt];

	if (c == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (rlo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	} else {
	  /* second or greater block of 16 */
#ifdef HAVE_AVX2
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[c-1][rlo-1],LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(matrix[c-1][rlo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	}

	debug15(print_vector_16(E_mask,rlo,c,"E_mask"));
	E_r_gap = _MM_MIN_EPI16(E_r_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E_r_gap"));
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_r, gap_open);
	dir_horiz = _MM_CMPGT_EPI16(E_r_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(T1,rlo,c,"T1"));
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_Egap"));

	E_r_gap = _MM_MAX_EPI16(E_r_gap, T1); /* Compare H + open with vert */
	E_r_gap = _MM_ADDS_EPI16(E_r_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_r_gap = _MM_MIN_EPI16(E_r_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_r_gap,rlo,c,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_r,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_r,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_r = _mm256_slli_si256(H_nogap_r,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_r,LAST_SHORT_SHIFT);
	H_nogap_r = _mm_slli_si128(H_nogap_r,ONE_SHORT);
#endif
	H_nogap_r = _MM_OR_SI(H_nogap_r, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores, allowing for alternate genomic nt */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm256_load_si256((__m256i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm256_adds_epi16(H_nogap_r, _mm256_max_epi16(pairscores_std,pairscores_alt));
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_std_ptr[rlo]));
	pairscores_alt = _mm_load_si128((__m128i *) &(pairscores_alt_ptr[rlo]));
	H_nogap_r = _mm_adds_epi16(H_nogap_r, _mm_max_epi16(pairscores_std,pairscores_alt));
#endif
	_mm_clflush(&H_nogap_r); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H"));

	dir_horiz = _MM_CMPGT_EPI16(E_r_gap,H_nogap_r); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[c][rlo]),dir_horiz);
#endif
	debug15(print_vector_16(dir_horiz,rlo,c,"dir_nogap"));

	H_nogap_r = _MM_MAX_EPI16(H_nogap_r, E_r_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_r,rlo,c,"H_nogap_r store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[rlo]), H_nogap_r);
#else
	_mm_store_si128((__m128i *) &(score_column[rlo]), H_nogap_r);
#endif

	/* Fix gaps along diagonal to prevent going into lower triangle, which can happen with ties between E and H */
	if (rhigh >= c) {
	  (*directions_Egap)[c][c] = DIAG;
	  (*directions_nogap)[c][c] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi16(E_mask,7);
	E_mask = _mm256_slli_si256(E_mask,ONE_SHORT);
	E_mask = _mm256_insert_epi16(E_mask,save,8);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
#endif
      }
    }
  }

#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = HORIZ;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16_upper\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		    revp,uband,/*upperp*/true);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
			rlength,glength,rsequence,gsequence,gsequence_alt,
			revp,uband,/*upperp*/true);
#endif

#ifdef CHECK1
  /* Check for row 0 directions */
  for (c = 1; c <= uband && c <= glength; c++) {
    assert((*directions_Egap)[c][0] != DIAG);
    assert((*directions_nogap)[c][0] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_16_upper_nonavx2(&directions_nogap_std,&directions_Egap_std,
					     this,rsequence,gsequence,gsequence_alt,
					     rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
					     open,extend,uband,jump_late_p,revp);
#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,
				goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,/*lband*/0,uband,jump_late_p,revp,/*saturation*/NEG_INFINITY_16,
				/*upperp*/true,/*lowerp*/false);
#endif

#ifdef DEBUG2
  printf("Banded\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
                    revp,uband,/*upperp*/true);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
                         rlength,glength,rsequence,gsequence,gsequence_alt,revp,uband,/*upperp*/true);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix16_compare_upper(matrix,matrix_std,rlength,glength,uband,
				rsequence,gsequence,gsequence_alt,
				goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions16_compare_nogap_upper(*directions_nogap,directions_nogap_std,rlength,glength,uband);

  banded_directions16_compare_Egap_upper(*directions_Egap,directions_Egap_std,rlength,glength,uband);
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef DEBUG_AVX2
/* Designed for computation below the diagonal, so no F loop or bottom masking needed */
/* Operates by rows */
Score16_T **
Dynprog_simd_16_lower_nonavx2 (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
			       T this, char *rsequence, char *gsequence, char *gsequence_alt,
			       int rlength, int glength,
			       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
			       Mismatchtype_T mismatchtype, int open, int extend,
			       int lband, bool jump_late_p, bool revp) {
  Score16_T **matrix, *score_column;
  __m128i pairscores_std;
  __m128i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, E_infinity, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_vert;
  int glength_ceil, r, c;
  int clo, chigh;
  int na1, na2, na2_alt;
  Score16_T *pairscores[5], *pairscores_ptr;
  Pairdistance_T **pairdistance_array_type, score1, score2;

  debug2(printf("Dynprog_simd_16_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  glength_ceil = (int) ((glength + SIMD_NSHORTS_NONAVX2)/SIMD_NSHORTS_NONAVX2) * SIMD_NSHORTS_NONAVX2;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_byrows (lower): "));
  debug(printf("Lengths are %d and %d, so band is %d on left\n",rlength,glength,lband));
  debug(printf("Genome length rounded up to %d\n",glength_ceil));

  matrix = aligned_score16_alloc(glength_ceil,rlength,
				 this->aligned_std.two.lower_matrix_ptrs,this->aligned_std.two.lower_matrix_space);
  *directions_nogap = aligned_directions16_alloc(glength_ceil,rlength,
						 this->aligned_std.two.lower_directions_ptrs_0,this->aligned_std.two.lower_directions_space_0);
  *directions_Egap = aligned_directions16_alloc(glength_ceil,rlength,
						this->aligned_std.two.lower_directions_ptrs_1,this->aligned_std.two.lower_directions_space_1);

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Egap)[r][0] = VERT;
    (*directions_nogap)[r][0] = VERT;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= SIMD_NSHORTS_NONAVX2 && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_nogap)[0][c] = HORIZ;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),16);
  pairscores[1] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),16);
  pairscores[2] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),16);
  pairscores[3] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),16);
  pairscores[4] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),16);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,glength_ceil*sizeof(Score16_T));
#endif

  c = 0; na2 = na2_alt = 'N';
  pairscores[0][c] = (Score16_T) pairdistance_array_type[(int) 'A'][na2];
  pairscores[1][c] = (Score16_T) pairdistance_array_type[(int) 'C'][na2];
  pairscores[2][c] = (Score16_T) pairdistance_array_type[(int) 'G'][na2];
  pairscores[3][c] = (Score16_T) pairdistance_array_type[(int) 'T'][na2];
  pairscores[4][c] = (Score16_T) pairdistance_array_type[(int) 'N'][na2];

  if (revp == false) {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[c-1];
      na2_alt = gsequence_alt[c-1];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score16_T) (score1 > score2) ? score1 : score2;
    }
  } else {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[1-c];
      na2_alt = gsequence_alt[1-c];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score16_T) (score1 > score2) ? score1 : score2;
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][c]),0,(glength_ceil-c)*sizeof(Score16_T));
#endif

  complement_dummy = _mm_set1_epi16(-1);

  gap_open = _mm_set1_epi16((Score16_T) open);
  gap_extend = _mm_set1_epi16((Score16_T) extend);

  E_infinity = _mm_set1_epi16(POS_INFINITY_16);
  if (jump_late_p) {
    for (clo = 0; clo <= glength; clo += SIMD_NSHORTS_NONAVX2) {
      if ((chigh = clo + SIMD_NSHORTS_NONAVX2 - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E >= H.  To fill in first row of each
	 column block with non-diags, make E == H. */
      E_mask = _mm_set1_epi16(1);

      E_c_gap = _mm_set1_epi16(NEG_INFINITY_16);
      H_nogap_c = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[rsequence[1-r]] : nt_to_int_array[rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	}

	debug15(print_vector_16(E_mask,clo,r,"E_mask"));
	E_c_gap = _mm_min_epi16(E_c_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_c, gap_open);
	dir_vert = _mm_cmplt_epi16(E_c_gap,T1); /* E < H */
	dir_vert = _mm_andnot_si128(dir_vert,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
	debug15(print_vector_16(T1,clo,r,"T1"));
	debug15(print_vector_16(dir_vert,clo,r,"dir_Egap"));

	E_c_gap = _mm_max_epi16(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _mm_adds_epi16(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_min_epi16(E_c_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_c,LAST_SHORT_NONAVX2);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_SHORT);
	H_nogap_c = _mm_or_si128(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence. */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm_adds_epi16(H_nogap_c, pairscores_std);
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_c,clo,r,"H"));

	dir_vert = _mm_cmplt_epi16(E_c_gap,H_nogap_c); /* E < H */
	dir_vert = _mm_andnot_si128(dir_vert,complement_dummy);	/* E >= H, for jump late */
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
	debug15(print_vector_16(dir_vert,clo,r,"dir_nogap"));

	H_nogap_c = _mm_max_epi16(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c store"));
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
      }
    }

  } else {
    /* jump early */
    for (clo = 0; clo <= glength; clo += SIMD_NSHORTS_NONAVX2) {
      if ((chigh = clo + SIMD_NSHORTS_NONAVX2 - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E > H.  To fill in first row of each
	 column block with non-diags, make E > H. */
      E_mask = _mm_set1_epi16(1);

      E_c_gap = _mm_set1_epi16(NEG_INFINITY_16+1);
      H_nogap_c = _mm_set1_epi16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[rsequence[1-r]] : nt_to_int_array[rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _mm_set1_epi16(0);
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _mm_set1_epi16(0);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
#endif
	} else {
	  /* second or greater block of 16 */
	  X_prev_nogap = _mm_set1_epi16(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_NONAVX2);
	}

	debug15(print_vector_16(E_mask,clo,r,"E_mask"));
	E_c_gap = _mm_min_epi16(E_c_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _mm_adds_epi16(H_nogap_c, gap_open);
	dir_vert = _mm_cmpgt_epi16(E_c_gap,T1); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
	debug15(print_vector_16(T1,clo,r,"T1"));
	debug15(print_vector_16(dir_vert,clo,r,"dir_Egap"));

	E_c_gap = _mm_max_epi16(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _mm_adds_epi16(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _mm_min_epi16(E_c_gap,_mm_add_epi16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E"));


	/* NOGAP */
	T1 = _mm_srli_si128(H_nogap_c,LAST_SHORT_NONAVX2);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_SHORT);
	H_nogap_c = _mm_or_si128(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm_adds_epi16(H_nogap_c, pairscores_std);
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_c,clo,r,"H"));

	dir_vert = _mm_cmpgt_epi16(E_c_gap,H_nogap_c); /* E > H, for jump early */
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
	debug15(print_vector_16(dir_vert,clo,r,"dir_nogap"));

	H_nogap_c = _mm_max_epi16(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c store"));
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = VERT;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16_lower\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		    revp,lband,/*upperp*/false);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
			rlength,glength,rsequence,gsequence,gsequence_alt,
			revp,lband,/*upperp*/false);
#endif

#ifdef CHECK1
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Egap)[r][0] != DIAG);
    assert((*directions_nogap)[r][0] != DIAG);
  }
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef HAVE_SSE2
/* Designed for computation below the diagonal, so no F loop or bottom masking needed */
/* Operates by rows */
Score16_T **
Dynprog_simd_16_lower (Direction16_T ***directions_nogap, Direction16_T ***directions_Egap,
		       T this, char *rsequence, char *gsequence, char *gsequence_alt,
		       int rlength, int glength,
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
		       int goffset, Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp,
#endif
		       Mismatchtype_T mismatchtype, int open, int extend,
		       int lband, bool jump_late_p, bool revp) {
  Score16_T **matrix, *score_column;
#ifdef HAVE_AVX2
  __m256i pairscores_std;
  __m256i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, E_infinity, T1;
  __m256i gap_open, gap_extend, complement_dummy;
  __m256i dir_vert;
  Score16_T save;
#else
  __m128i pairscores_std;
  __m128i H_nogap_c, X_prev_nogap, E_c_gap, E_mask, E_infinity, T1;
  __m128i gap_open, gap_extend, complement_dummy;
  __m128i dir_vert;
#endif
  int glength_ceil, r, c;
  int clo, chigh;
  int na1, na2, na2_alt;
  Score16_T *pairscores[5], *pairscores_ptr;
  Pairdistance_T **pairdistance_array_type, score1, score2;

#ifdef DEBUG_AVX2
  Score16_T **matrix_std;
  Direction16_T **directions_nogap_std, **directions_Egap_std;
  char na2_single;
#elif defined(DEBUG_SIMD)
  Score32_T **matrix_std;
  Direction32_T **directions_nogap_std, **directions_Egap_std, **directions_Fgap_std;
  char na2_single;
#endif


  debug2(printf("Dynprog_simd_16_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));
  debug15(printf("Dynprog_simd_16_lower.  jump_late_p %d, open %d, extend %d\n",jump_late_p,open,extend));

  glength_ceil = (int) ((glength + SIMD_NSHORTS)/SIMD_NSHORTS) * SIMD_NSHORTS;
  pairdistance_array_type = pairdistance_array[mismatchtype];
  
  debug(printf("compute_scores_simd_16_byrows (lower): "));
  debug(printf("Lengths are %d and %d, so band is %d on left\n",rlength,glength,lband));
  debug(printf("Genome length rounded up to %d\n",glength_ceil));

  matrix = aligned_score16_alloc(glength_ceil,rlength,
				 this->aligned.two.lower_matrix_ptrs,this->aligned.two.lower_matrix_space);
  *directions_nogap = aligned_directions16_alloc(glength_ceil,rlength,
						 this->aligned.two.lower_directions_ptrs_0,this->aligned.two.lower_directions_space_0);
  *directions_Egap = aligned_directions16_alloc(glength_ceil,rlength,
						this->aligned.two.lower_directions_ptrs_1,this->aligned.two.lower_directions_space_1);

#if 0
  /* Column 0 initialization */
  /* penalty = open; */
  for (r = 1; r <= lband && r <= rlength; r++) {
    /* penalty += extend; */
    (*directions_Egap)[r][0] = VERT;
    (*directions_nogap)[r][0] = VERT;
  }
#endif
#if 0
  /* Already initialized to DIAG.  Actually, no longer initializing directions_Egap */
  (*directions_Egap)[1][0] = DIAG; /* previously used STOP */
  (*directions_nogap)[0][0] = DIAG; /* previously used STOP */
#endif
#if 0
  /* Row 0 initialization */
  /* penalty = open; */
  for (c = 1; c <= SIMD_NSHORTS && c <= glength; c++) {
    /* penalty += extend; */
    (*directions_nogap)[0][c] = HORIZ;
  }
#endif


  /* Load pairscores.  Store match - mismatch */
  pairscores[0] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[1] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[2] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[3] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),ALIGN_SIZE);
  pairscores[4] = (Score16_T *) _mm_malloc(glength_ceil * sizeof(Score16_T),ALIGN_SIZE);

#if 0
  /* Should not be necessary */
  memset((void *) pairscores[0],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[1],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[2],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[3],0,glength_ceil*sizeof(Score16_T));
  memset((void *) pairscores[4],0,glength_ceil*sizeof(Score16_T));
#endif

  c = 0; na2 = na2_alt = 'N';
  pairscores[0][c] = (Score16_T) pairdistance_array_type[(int) 'A'][na2];
  pairscores[1][c] = (Score16_T) pairdistance_array_type[(int) 'C'][na2];
  pairscores[2][c] = (Score16_T) pairdistance_array_type[(int) 'G'][na2];
  pairscores[3][c] = (Score16_T) pairdistance_array_type[(int) 'T'][na2];
  pairscores[4][c] = (Score16_T) pairdistance_array_type[(int) 'N'][na2];

  if (revp == false) {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[c-1];
      na2_alt = gsequence_alt[c-1];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score16_T) (score1 > score2) ? score1 : score2;
    }
  } else {
    for (c = 1; c <= glength; c++) {
      na2 = gsequence[1-c];
      na2_alt = gsequence_alt[1-c];
      /* Take max here */
      score1 = pairdistance_array_type[(int) 'A'][na2];
      score2 = pairdistance_array_type[(int) 'A'][na2_alt];
      pairscores[0][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'C'][na2];
      score2 = pairdistance_array_type[(int) 'C'][na2_alt];
      pairscores[1][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'G'][na2];
      score2 = pairdistance_array_type[(int) 'G'][na2_alt];
      pairscores[2][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'T'][na2];
      score2 = pairdistance_array_type[(int) 'T'][na2_alt];
      pairscores[3][c] = (Score16_T) (score1 > score2) ? score1 : score2;

      score1 = pairdistance_array_type[(int) 'N'][na2];
      score2 = pairdistance_array_type[(int) 'N'][na2_alt];
      pairscores[4][c] = (Score16_T) (score1 > score2) ? score1 : score2;
    }
  }

#if 0
  /* Should not be necessary */
  memset((void *) &(pairscores[0][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[1][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[2][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[3][c]),0,(glength_ceil-c)*sizeof(Score16_T));
  memset((void *) &(pairscores[4][c]),0,(glength_ceil-c)*sizeof(Score16_T));
#endif

  complement_dummy = _MM_SET1_EPI16(-1);

  gap_open = _MM_SET1_EPI16((Score16_T) open);
  gap_extend = _MM_SET1_EPI16((Score16_T) extend);

  E_infinity = _MM_SET1_EPI16(POS_INFINITY_16);
  if (jump_late_p) {
    for (clo = 0; clo <= glength; clo += SIMD_NSHORTS) {
      if ((chigh = clo + SIMD_NSHORTS - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E >= H.  To fill in first row of each
	 column block with non-diags, make E == H. */
      E_mask = _MM_SET1_EPI16(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_c_gap = _MM_SET1_EPI16(NEG_INFINITY_16);
      H_nogap_c = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[(int) rsequence[1-r]] : nt_to_int_array[(int) rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	} else {
	  /* second or greater block of 16 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[r-1][clo-1],LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	}

	debug15(print_vector_16(E_mask,clo,r,"E_mask"));
	E_c_gap = _MM_MIN_EPI16(E_c_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_c, gap_open);
	dir_vert = _MM_CMPLT_EPI16(E_c_gap,T1); /* E < H */
	dir_vert = _MM_ANDNOT_SI(dir_vert,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_16(T1,clo,r,"T1"));
	debug15(print_vector_16(dir_vert,clo,r,"dir_Egap"));

	E_c_gap = _MM_MAX_EPI16(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _MM_ADDS_EPI16(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _MM_MIN_EPI16(E_c_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_c,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_c,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_c = _mm256_slli_si256(H_nogap_c,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_c,LAST_SHORT_SHIFT);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_SHORT);
#endif
	H_nogap_c = _MM_OR_SI(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence. */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm256_adds_epi16(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm_adds_epi16(H_nogap_c, pairscores_std);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_c,clo,r,"H"));

	dir_vert = _MM_CMPLT_EPI16(E_c_gap,H_nogap_c); /* E < H */
	dir_vert = _MM_ANDNOT_SI(dir_vert,complement_dummy);	/* E >= H, for jump late */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_16(dir_vert,clo,r,"dir_nogap"));

	H_nogap_c = _MM_MAX_EPI16(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[clo]), H_nogap_c);
#else
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);
#endif


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi16(E_mask,7);
	E_mask = _mm256_slli_si256(E_mask,ONE_SHORT);
	E_mask = _mm256_insert_epi16(E_mask,save,8);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
#endif
      }
    }

  } else {
    /* jump early */
    for (clo = 0; clo <= glength; clo += SIMD_NSHORTS) {
      if ((chigh = clo + SIMD_NSHORTS - 1) > glength) {
	chigh = glength;
      }

      /* dir_vert tests if E > H.  To fill in first row of each
	 column block with non-diags, make E > H. */
      E_mask = _MM_SET1_EPI16(1);

      /* Holds for all INITIAL_GAP_PENALTY */
      E_c_gap = _MM_SET1_EPI16(NEG_INFINITY_16+1);
      H_nogap_c = _MM_SET1_EPI16(NEG_INFINITY_16-open); /* Compensate for T1 = H + open */

      for (r = clo; r <= chigh + lband && r <= rlength; r++) {
	score_column = matrix[r];

	if (r == 0) {
	  na1 = 4; /* 'N' */
	} else {
	  na1 = revp ? nt_to_int_array[(int) rsequence[1-r]] : nt_to_int_array[(int) rsequence[r-1]];
	}
	pairscores_ptr = pairscores[na1];

	if (r == 0) {
	  X_prev_nogap = _MM_SETZERO_SI();
	} else if (clo == 0) {
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),NEG_INFINITY_16,LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(NEG_INFINITY_16); /* works if we start outside the rlo bounds */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	} else {
	  /* second or greater block of 16 */
#ifdef ZERO_INITIAL_GAP_PENALTY
	  X_prev_nogap = _MM_SETZERO_SI();
#elif defined(HAVE_AVX2)
	  X_prev_nogap = _mm256_insert_epi16(_mm256_setzero_si256(),matrix[r-1][clo-1],LAST_SHORT_INSERT);
#else
	  X_prev_nogap = _mm_set1_epi16(matrix[r-1][clo-1]); /* get H from previous block and previous column */
	  X_prev_nogap = _mm_srli_si128(X_prev_nogap,LAST_SHORT_SHIFT);
#endif
	}

	debug15(print_vector_16(E_mask,clo,r,"E_mask"));
	E_c_gap = _MM_MIN_EPI16(E_c_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E_c_gap"));
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c load"));

	/* EGAP */
	T1 = _MM_ADDS_EPI16(H_nogap_c, gap_open);
	dir_vert = _MM_CMPGT_EPI16(E_c_gap,T1); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_Egap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_Egap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_16(T1,clo,r,"T1"));
	debug15(print_vector_16(dir_vert,clo,r,"dir_Egap"));

	E_c_gap = _MM_MAX_EPI16(E_c_gap, T1); /* Compare H + open with vert */
	E_c_gap = _MM_ADDS_EPI16(E_c_gap, gap_extend); /* Compute scores for Egap (vert + open) */
	E_c_gap = _MM_MIN_EPI16(E_c_gap,_MM_ADD_EPI16(E_mask,E_infinity));
	debug15(print_vector_16(E_c_gap,clo,r,"E"));


	/* NOGAP */
#ifdef HAVE_AVX2
	T1 = _mm256_insert_epi16(_mm256_setzero_si256(),_mm256_extract_epi16(H_nogap_c,SIMD_NSHORTS-1),LAST_SHORT_INSERT);
	X_prev_nogap = _mm256_insert_epi16(X_prev_nogap,_mm256_extract_epi16(H_nogap_c,MID_SHORT_INSERT-1),MID_SHORT_INSERT);
	H_nogap_c = _mm256_slli_si256(H_nogap_c,ONE_SHORT);
#else
	T1 = _mm_srli_si128(H_nogap_c,LAST_SHORT_SHIFT);
	H_nogap_c = _mm_slli_si128(H_nogap_c,ONE_SHORT);
#endif
	H_nogap_c = _MM_OR_SI(H_nogap_c, X_prev_nogap);
	X_prev_nogap = T1;

	/* Add pairscores.  No alternate chars for query sequence */
#ifdef HAVE_AVX2
	pairscores_std = _mm256_load_si256((__m256i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm256_adds_epi16(H_nogap_c, pairscores_std);
#else
	pairscores_std = _mm_load_si128((__m128i *) &(pairscores_ptr[clo]));
	H_nogap_c = _mm_adds_epi16(H_nogap_c, pairscores_std);
#endif
	_mm_clflush(&H_nogap_c); /* Needed for opencc -O3 on AMD */
	debug15(print_vector_16(H_nogap_c,clo,r,"H"));

	dir_vert = _MM_CMPGT_EPI16(E_c_gap,H_nogap_c); /* E > H, for jump early */
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &((*directions_nogap)[r][clo]),dir_vert);
#else
	_mm_store_si128((__m128i *) &((*directions_nogap)[r][clo]),dir_vert);
#endif
	debug15(print_vector_16(dir_vert,clo,r,"dir_nogap"));

	H_nogap_c = _MM_MAX_EPI16(H_nogap_c, E_c_gap); /* Compare H + pairscores with horiz + extend */
	debug15(print_vector_16(H_nogap_c,clo,r,"H_nogap_c store"));
#ifdef HAVE_AVX2
	_mm256_store_si256((__m256i *) &(score_column[clo]), H_nogap_c);
#else
	_mm_store_si128((__m128i *) &(score_column[clo]), H_nogap_c);
#endif


	/* Fix gaps along diagonal to prevent going into upper triangle, which can happen with ties between E and H */
	if (chigh >= r) {
	  (*directions_Egap)[r][r] = DIAG;
	  (*directions_nogap)[r][r] = DIAG;
	}

	/* No need for F loop here */
#ifdef HAVE_AVX2
	save = _mm256_extract_epi16(E_mask,7);
	E_mask = _mm256_slli_si256(E_mask,ONE_SHORT);
	E_mask = _mm256_insert_epi16(E_mask,save,8);
#else
	E_mask = _mm_slli_si128(E_mask,ONE_SHORT);
#endif
      }
    }
  }


#ifdef CHECK1
  /* Row 0 and column 0 directions fail anyway due to saturation */
  /* Handle (0,1) and (1,0) directions, otherwise DIAG */
  (*directions_Egap)[1][0] = VERT;
#endif

#ifdef DEBUG2
  printf("SIMD: Dynprog_simd_16_lower\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		    revp,lband,/*upperp*/false);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
			rlength,glength,rsequence,gsequence,gsequence_alt,
			revp,lband,/*upperp*/false);
#endif

#ifdef CHECK1
  /* Check for column 0 directions */
  for (r = 1; r <= lband && r <= rlength; r++) {
    assert((*directions_Egap)[r][0] != DIAG);
    assert((*directions_nogap)[r][0] != DIAG);
  }
#endif

#ifdef DEBUG_AVX2
  matrix_std = Dynprog_simd_16_lower_nonavx2(&directions_nogap_std,&directions_Egap_std,
					     this,rsequence,gsequence,gsequence_alt,
					     rlength,glength,goffset,chroffset,chrhigh,watsonp,mismatchtype,
					     open,extend,lband,jump_late_p,revp);
#elif defined(DEBUG_SIMD)
  matrix_std = Dynprog_standard(&directions_nogap_std,&directions_Egap_std,&directions_Fgap_std,
				this,rsequence,/*gsequence (NULL for debugging)*/NULL,/*gsequence_alt*/NULL,
				rlength,glength,
				goffset,chroffset,chrhigh,watsonp,mismatchtype,
				open,extend,lband,/*uband*/0,jump_late_p,revp,/*saturation*/NEG_INFINITY_16,
				/*upperp*/false,/*lowerp*/true);
#endif

#ifdef DEBUG2
  printf("Banded\n");
  Matrix16_print_ud(matrix,rlength,glength,rsequence,gsequence,gsequence_alt,
		    revp,lband,/*upperp*/false);
  Directions16_print_ud(*directions_nogap,*directions_Egap,
			rlength,glength,rsequence,gsequence,gsequence_alt,
			revp,lband,/*upperp*/false);
#endif
  
#if defined(DEBUG_AVX2) || defined(DEBUG_SIMD)
  banded_matrix16_compare_lower(matrix,matrix_std,rlength,glength,lband,
				rsequence,gsequence,gsequence_alt,
				goffset,chroffset,chrhigh,watsonp,revp);

  banded_directions16_compare_nogap_lower(*directions_nogap,directions_nogap_std,rlength,glength,lband);

  banded_directions16_compare_Egap_lower(*directions_Egap,directions_Egap_std,rlength,glength,lband);
#endif

  _mm_free(pairscores[4]);
  _mm_free(pairscores[3]);
  _mm_free(pairscores[2]);
  _mm_free(pairscores[1]);
  _mm_free(pairscores[0]);

  return matrix;
}
#endif


#ifdef DEBUG17
static char complCode[128] = COMPLEMENT_LC;

static char
get_genomic_nt (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

#if 0
  /* If the read has a deletion, then we will extend beyond 0 or genomiclength, so do not restrict. */
  if (genomicpos < 0) {
    return '*';

  } else if (genomicpos >= genomiclength) {
    return '*';

  }
#endif

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char(genome,pos)));
      return Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      debug8(printf("At %u, genomicnt is %c\n",
		    genomicpos,Genome_get_char_blocks(pos)));
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    if ((pos = chrhigh - genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';

    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';

#if 0
    } else if (genome) {
      /* Not necessary, because Genome_get_char_blocks should work */
      c2 = Genome_get_char(genome,pos);
#endif

    } else {
      /* GMAP with user-supplied genomic segment */
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    debug8(printf("At %u, genomicnt is %c\n",genomicpos,complCode[(int) c2]));
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
List_T
Dynprog_traceback_8 (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		     Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		     int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
		     int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		     Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction8_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_8 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) == HORIZ) {
      dist = 1;
      while (c > 0 && directions_Egap[c--][r] != DIAG) {
	dist++;
      }
#if 0
      if (c == 0) {
	/* Directions in column 0 can sometimes be DIAG */
	dir = VERT;
      } else {
	printf("| ");		/* For Fgap */
	dir = directions_nogap[c][r];
      }
#endif

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,/*genomesequence*/NULL,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      watsonp,dynprogindex);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else if (dir == VERT) {
      dist = 1;
      while (r > 0 && directions_Fgap[c][r--] != DIAG) {
	dist++;
      }
#if 0
      if (r == 0) {
	/* Directions in row 0 can sometimes be DIAG */
	dir = HORIZ;
      } else {
	dir = directions_nogap[c][r];
      }
#endif

      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,r+dist,c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else if (dir == DIAG) {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));
	
      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
      
      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);
      
      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;

    } else {
      fprintf(stderr,"Bad dir at r %d, c %d\n",r,c);
      abort();
    }
  }

  if (r == 0 && c == 0) {
    /* Finished with a diagonal step */

  } else if (c == 0) {
    dist = r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,r,/*c*/0+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));

  } else {
    assert(r == 0);
    dist = c;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/*r*/0+LAZY_INDEL,c,dist,/*genomesequence*/NULL,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    watsonp,dynprogindex);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
List_T
Dynprog_traceback_8_upper (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			   Direction8_T **directions_nogap, Direction8_T **directions_Egap,
			   int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
			   int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			   Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction8_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_8_upper at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) != DIAG) {
      /* Must be HORIZ */
      dist = 1;
      /* Should not need to check for c > 0 if the main diagonal is populated with DIAG */
      while (/* c > 0 && */ directions_Egap[c--][r] != DIAG) {
	dist++;
      }
      /* assert(c != 0); */

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,/*genomesequence*/NULL,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      watsonp,dynprogindex);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));
	
      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);
	
      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;
    }
  }

  assert(r == 0);
  if (/* r == 0 && */ c == 0) {
    /* Finished with a diagonal step */

  } else {
    assert(c != 0);
    assert(r == 0);
    dist = c;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/*r*/0+LAZY_INDEL,c,dist,/*genomesequence*/NULL,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    watsonp,dynprogindex);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}

List_T
Dynprog_traceback_8_lower (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			   Direction8_T **directions_nogap, Direction8_T **directions_Egap,
			   int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
			   int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  int querycoord, genomecoord;
  Direction8_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_8_lower at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[r][c]) != DIAG) {
      /* Must be VERT */
      dist = 1;
      /* Should not need to check for r > 0 if the main diagonal is populated with DIAG */
      while (/* r > 0 && */ directions_Egap[r--][c] != DIAG) {
	dist++;
      }
      /* assert(r != 0); */

      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,r+dist,c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));
	
      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;
    }
  }

  assert(c == 0);
  if (r == 0 /* && c == 0 */) {
    /* Finished with a diagonal step */

  } else {
    assert(r != 0);
    assert(c == 0);
    dist = r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,r,/*c*/0+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  return pairs;
}


List_T
Dynprog_traceback_16 (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
		      Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		      int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
		      int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
		      Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction16_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_16 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) == HORIZ) {
      dist = 1;
      while (c > 0 && directions_Egap[c--][r] != DIAG) {
	dist++;
      }
#if 0
      if (c == 0) {
	/* Directions in column 0 can sometimes be DIAG */
	dir = VERT;
      } else {
	dir = directions_nogap[c][r];
      }
#endif

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,/*genomesequence*/NULL,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      watsonp,dynprogindex);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));
      
    } else if (dir == VERT) {
      dist = 1;
      while (r > 0 && directions_Fgap[c][r--] != DIAG) {
	dist++;
      }
#if 0
      if (r == 0) {
	/* Directions in row 0 can sometimes be DIAG */
	dir = HORIZ;
      } else {
	dir = directions_nogap[c][r];
      }
#endif

      debug(printf("V%d: ",dist));
      debug(printf("New dir at %d,%d is %d\n",c,r,dir));
      pairs = Pairpool_add_queryskip(pairs,r+dist,c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else if (dir == DIAG) {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;

    } else {
      fprintf(stderr,"Bad dir at r %d, c %d\n",r,c);
      abort();
    }
  }

  if (r == 0 && c == 0) {
    /* Finished with a diagonal step */

  } else if (c == 0) {
    dist = r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,r,/*c*/0+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));

  } else {
    assert(r == 0);
    dist = c;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/*r*/0+LAZY_INDEL,c,dist,/*genomesequence*/NULL,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    watsonp,dynprogindex);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}


List_T
Dynprog_traceback_16_upper (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			    Direction16_T **directions_nogap, Direction16_T **directions_Egap,
			    int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
			    int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  bool add_dashes_p;
  int querycoord, genomecoord;
  Direction16_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_16_upper at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) != DIAG) {
      /* Must be HORIZ */
      dist = 1;
      /* Should not need to check for c > 0 if the main diagonal is populated with DIAG */
      while (/* c > 0 && */ directions_Egap[c--][r] != DIAG) {
	dist++;
      }
      /* assert(c != 0); */

      debug(printf("H%d: ",dist));
      pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,r,c+dist,dist,/*genomesequence*/NULL,
				      queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				      watsonp,dynprogindex);
      if (add_dashes_p == true) {
	*nopens += 1;
	*nindels += dist;
      }
      debug(printf("\n"));

    } else {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;
    }
  }

  assert(r == 0);
  if (/* r == 0 && */ c == 0) {
    /* Finished with a diagonal step */

  } else {
    assert(c != 0);
    assert(r == 0);
    dist = c;
    debug(printf("H%d: ",dist));
    pairs = Pairpool_add_genomeskip(&add_dashes_p,pairs,/*r*/0+LAZY_INDEL,c,dist,/*genomesequence*/NULL,
				    queryoffset,genomeoffset,pairpool,revp,chroffset,chrhigh,
				    watsonp,dynprogindex);
    if (add_dashes_p == true) {
      *nopens += 1;
      *nindels += dist;
    }
    debug(printf("\n"));
  }

  return pairs;
}

List_T
Dynprog_traceback_16_lower (List_T pairs, int *nmatches, int *nmismatches, int *nopens, int *nindels,
			    Direction16_T **directions_nogap, Direction16_T **directions_Egap,
			    int r, int c, char *rsequence, char *rsequenceuc, char *gsequence, char *gsequence_alt,
			    int queryoffset, int genomeoffset, Pairpool_T pairpool, bool revp, int dynprogindex) {
  char c1, c1_uc, c2, c2_alt;
  int dist;
  int querycoord, genomecoord;
  Direction16_T dir;
#ifdef DEBUG17
  char c2_single;
#endif

  debug(printf("Starting traceback_16_lower at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  while (r > 0 && c > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[r][c]) != DIAG) {
      /* Must be VERT */
      dist = 1;
      /* Should not need to check for r > 0 if the main diagonal is populated with DIAG */
      while (/* r > 0 && */ directions_Egap[r--][c] != DIAG) {
	dist++;
      }
      /* assert(r != 0); */

      debug(printf("V%d: ",dist));
      pairs = Pairpool_add_queryskip(pairs,r+dist,c,dist,rsequence,
				     queryoffset,genomeoffset,pairpool,revp,
				     dynprogindex);
      *nopens += 1;
      *nindels += dist;
      debug(printf("\n"));

    } else {
      querycoord = r-1;
      genomecoord = c-1;
      if (revp == true) {
	querycoord = -querycoord;
	genomecoord = -genomecoord;
      }

      c1 = rsequence[querycoord];
      c1_uc = rsequenceuc[querycoord];
      c2 = gsequence[genomecoord];
      c2_alt = gsequence_alt[genomecoord];
#ifdef DEBUG17
      c2_single = get_genomic_nt(&c2_alt,genomeoffset+genomecoord,chroffset,chrhigh,watsonp);
      if (c2 != c2_single) {
	abort();
      }
#endif

#ifdef EXTRACT_GENOMICSEG
      assert(c2 == genomesequence[genomecoord]);
#endif

      if (c2 == '*') {
	/* Don't push pairs past end of chromosome */
	debug(printf("Don't push pairs past end of chromosome: genomeoffset %u, genomecoord %u, chroffset %u, chrhigh %u, watsonp %d\n",
		     genomeoffset,genomecoord,chroffset,chrhigh,watsonp));

      } else if (/*querysequenceuc[querycoord]*/c1_uc == c2 || c1_uc == c2_alt) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - match\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,DYNPROG_MATCH_COMP,c2,c2_alt,dynprogindex);

      } else if (consistent_array[(int) c1_uc][(int) c2] == true || consistent_array[(int) c1_uc][(int) c2_alt] == true) {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - ambiguous\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,AMBIGUOUS_COMP,c2,c2_alt,dynprogindex);

      } else {
	debug(printf("Pushing %d,%d [%d,%d] (%c,%c) - mismatch\n",
		     r,c,queryoffset+querycoord,genomeoffset+genomecoord,c1_uc,c2));
	*nmismatches += 1;
	pairs = Pairpool_push(pairs,pairpool,queryoffset+querycoord,genomeoffset+genomecoord,
			      c1,MISMATCH_COMP,c2,c2_alt,dynprogindex);
      }

      r--; c--;
    }
  }

  assert(c == 0);
  if (r == 0 /* && c == 0 */) {
    /* Finished with a diagonal step */

  } else {
    assert(r != 0);
    assert(c == 0);
    dist = r;
    debug(printf("V%d: ",dist));
    pairs = Pairpool_add_queryskip(pairs,r,/*c*/0+LAZY_INDEL,dist,rsequence,
				   queryoffset,genomeoffset,pairpool,revp,
				   dynprogindex);
    *nopens += 1;
    *nindels += dist;
    debug(printf("\n"));
  }

  return pairs;
}
#endif


