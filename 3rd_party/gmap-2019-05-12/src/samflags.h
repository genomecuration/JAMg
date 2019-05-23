/* $Id: samflags.h 214301 2018-03-19 23:36:38Z twu $ */
#ifndef SAMFLAGS_INCLUDED
#define SAMFLAGS_INCLUDED

#define PAIRED_READ        0x0001 /* 1 */
#define PAIRED_MAPPING     0x0002 /* 2 */
#define QUERY_UNMAPPED     0x0004 /* 4 */
#define MATE_UNMAPPED      0x0008 /* 8 */
#define QUERY_MINUSP       0x0010 /* 16 */
#define MATE_MINUSP        0x0020 /* 32 */
#define FIRST_READ_P       0x0040 /* 64 */
#define SECOND_READ_P      0x0080 /* 128 */
#define NOT_PRIMARY        0x0100 /* 256 */
#define BAD_READ_QUALITY   0x0200 /* 512 */
#define DUPLICATE_READ     0x0400 /* 1024 */
#define SUPPLEMENTARY      0x0800 /* 2048 */


/* 83 = first read, minus strand for paired */
/* 99 = first read, plus strand for paired */
/* 147 = second read, minus strand for paired */
/* 163 = second read, plus strand for paired */

/* For forcing a read to be primary */
#define SET_PRIMARY        0xFEFF /* do a logical-and (&) with this */


/* XO tag for output type */
#define ABBREV_NOMAPPING_1 "NM"
#define ABBREV_NOMAPPING_2 "NM"

#define ABBREV_UNPAIRED_UNIQ "UU"
#define ABBREV_UNPAIRED_TRANSLOC "UT"
#define ABBREV_UNPAIRED_MULT "UM"

#define ABBREV_UNPAIRED_CIRCULAR "UC"
#define ABBREV_UNPAIRED_MULT_XS "UX"


#define ABBREV_HALFMAPPING_UNIQ "HU"
#define ABBREV_HALFMAPPING_TRANSLOC "HT"
#define ABBREV_HALFMAPPING_MULT "HM"

#define ABBREV_PAIRED_UNIQ_INV "PI"
#define ABBREV_PAIRED_UNIQ_SCR "PS"
#define ABBREV_PAIRED_UNIQ_LONG "PL"
#define ABBREV_PAIRED_MULT "PM"

#define ABBREV_CONCORDANT_UNIQ "CU"
#define ABBREV_CONCORDANT_TRANSLOC "CT"
#define ABBREV_CONCORDANT_MULT "CM"

#define ABBREV_HALFMAPPING_CIRCULAR "HC"
#define ABBREV_PAIRED_UNIQ_CIRCULAR "PC"
#define ABBREV_CONCORDANT_CIRCULAR "CC"

#define ABBREV_HALFMAPPING_MULT_XS "HX"
#define ABBREV_PAIRED_MULT_XS "PX"
#define ABBREV_CONCORDANT_MULT_XS "CX"


typedef enum {OUTPUT_FILE,	/* 0: specified by the --output-file option */

	      OUTPUT_NONE,	/* 1: used when omit_concordant_uniq_p or omit_concordant_mult_p is set */

	      OUTPUT_NM,	/* 2: nomapping */

	      OUTPUT_UU,	/* 3: unpaired_uniq */
	      OUTPUT_UT,	/* 4: unpaired_transloc */
	      OUTPUT_UM,	/* 5: unpaired_mult (N_SPLIT_OUTPUTS_SINGLE_STD) */

	      OUTPUT_UC,	/* 6: unpaired_circular (N_SPLIT_OUTPUTS_SINGLE_TOCIRC) */
	      OUTPUT_UX,	/* 7: unpaired_mult_xs (N_SPLIT_OUTPUTS_SINGLE) */


	      OUTPUT_HU,	/* 8: halfmapping_uniq */
	      OUTPUT_HT,	/* 9: halfmapping_transloc */
	      OUTPUT_HM,	/* 10: halfmapping_mult */

	      OUTPUT_PI,	/* 11: paired_uniq_inv */
	      OUTPUT_PS,	/* 12: paired_uniq_scr */
	      OUTPUT_PL,	/* 13: paired_uniq_long */
	      OUTPUT_PM,	/* 14: paired_mult */
	      OUTPUT_CU,	/* 15: concordant_uniq */
	      OUTPUT_CT,	/* 16: concordant_transloc */
	      OUTPUT_CM,	/* 17: concordant_mult (N_SPLIT_OUTPUTS_STD) */

	      OUTPUT_HC,	/* 18: halfmapping_circular */
	      OUTPUT_PC,	/* 19: paired_uniq_circular */
	      OUTPUT_CC,	/* 20: concordant_circular (N_SPLIT_OUTPUTS_TOCIRC) */

	      OUTPUT_HX,	/* 21: halfmapping_mult_xs */
	      OUTPUT_PX,	/* 22: paired_mult_xs */
	      OUTPUT_CX}	/* 23: concordant_mult_xs (N_SPLIT_OUTPUTS) */
  SAM_split_output_type;

/* output 0 is stdout */
/* Code checks for split_output = 1 to N_SPLIT_OUTPUTS inclusive */
#define N_SPLIT_OUTPUTS_SINGLE_STD 5
#define N_SPLIT_OUTPUTS_SINGLE_TOCIRC 6
#define N_SPLIT_OUTPUTS_SINGLE 7

#define N_SPLIT_OUTPUTS_STD 17
#define N_SPLIT_OUTPUTS_TOCIRC 20
#define N_SPLIT_OUTPUTS 23




/* GSNAP outputs */
#if 0
  FILE *fp_nomapping;		/* NM */

  FILE *fp_unpaired_uniq;	 /* UU */
  FILE *fp_unpaired_transloc;	 /* UT */
  FILE *fp_unpaired_mult;	 /* UM */

  FILE *fp_unpaired_circular;	 /* UC */
  FILE *fp_unpaired_mult_xs_1;	 /* UX */
  FILE *fp_unpaired_mult_xs_2;	 /* UX */

  FILE *fp_halfmapping_uniq;	/* HU */
  FILE *fp_halfmapping_transloc; /* HT */
  FILE *fp_halfmapping_mult;	 /* HM */

  FILE *fp_paired_uniq_inv;	 /* PI */
  FILE *fp_paired_uniq_scr;	 /* PS */
  FILE *fp_paired_uniq_long;	 /* PL */
  FILE *fp_paired_mult;		 /* PM */

  FILE *fp_concordant_uniq;	 /* CU */
  FILE *fp_concordant_transloc;	 /* CT */
  FILE *fp_concordant_mult;	 /* CM */

  FILE *fp_halfmapping_circular; /* HC */
  FILE *fp_paired_uniq_circular; /* PC */
  FILE *fp_concordant_circular;	 /* CC */

  FILE *fp_halfmapping_mult_xs_1; /* HX */
  FILE *fp_halfmapping_mult_xs_2; /* HX */
  FILE *fp_paired_mult_xs_1;	 /* PX */
  FILE *fp_paired_mult_xs_2;	 /* PX */
  FILE *fp_concordant_mult_xs_1; /* CX */
  FILE *fp_concordant_mult_xs_2; /* CX */
#endif

/* GMAP outputs */
#if 0
  FILE *fp_nomapping;		/* NM */

  FILE *fp_uniq;		/* UU */
  FILE *fp_transloc;		/* UT */
  FILE *fp_mult;		/* UM */

  FILE *fp_circular;		/* UC */
  FILE *fp_mult_xs;		/* UX */
#endif


#endif

