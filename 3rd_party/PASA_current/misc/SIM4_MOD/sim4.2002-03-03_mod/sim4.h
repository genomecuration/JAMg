#ifndef SIM4_H
#define SIM4_H

/* $Id: sim4.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

#define  DEFAULT_NUM_I     3       
#define  DIST_CUTOFF       3
#define  DEFAULT_MIN_COV   (.8)

#define  MININT           (-99999)
#define  MIN_INTRON        30
#define  MAX_INTRON        20000
#define  MAX_GRINIT        500
#define  MAX_INTERNAL_GAP  50   /* 50 */
#define  LL                60

#define  DEFAULT_DRANGE    10
#define  DEFAULT_WEIGHT    100
#define  DEFAULT_RELINK_H  500
#define  DEFAULT_W         12
#define  DEFAULT_X         12    
#define  DEFAULT_K         16
#define  DEFAULT_C         12
#define  LINK              0
#define  RELINK            1
#define  P                 (.2)

#define  min(x,y)        ((x>y) ? (y):(x))
#define  max(x,y)        ((x<y) ? (y):(x))
#define  START_SIG       ((G_score >= abs(C_score)) ? "GT" : "CT") 
#define  END_SIG         ((G_score >= abs(C_score)) ? "AG" : "AC")

#define  MATCH           1
#define  MISMATCH       (-5)
#define  L               8

#define  DELETE          1
#define  INSERT          2
#define  SUBSTITUTE      3
#define  INTRON          4
#define  O_INTRON        5

#undef TRUE
#undef FALSE
enum { FALSE = 0, TRUE = 1};
enum { INIT = 0, PERM = 1, TEMP = 2};
enum { EST_GEN = 1, GEN_EST = 2 };
enum { FWD = 0, BWD = 1, BOTH = 2 };
enum { OK = 0, FREE_START = 1, FREE_END = 2, FREE_BOTH_ENDS = 3};

/* data structures */

/* used in select_path() */
typedef struct msp {
        int len, pos1, pos2;
        int score;
        int Score;
        int  prev;
        struct msp *next_msp;
        } *Msp_ptr;

typedef struct exon {
        int  from1, from2, to1, to2;
        int  min_diag, max_diag;
        int  match;
        char ori;
        int  length;
        int  flag;
        int  ematches;
        int  nmatches;
        int  edist;
        int  alen;
        Msp_ptr msps;
        struct exon *next_exon;
        } *Exon_ptr;

typedef struct intron {
        int from1, from2, to1, to2;
        int  length;
        char orientation;
        struct intron *next_intron;
} *Intron_ptr, Intron;

typedef struct exon   Exon;

typedef struct coordinates {
        int pos1;
        int pos2;
}  coords;

/* used only in the alignment stage */
typedef struct edit_script {
        char op_type;   /* SUB, INS, or DEL */
        int num;        /* Number of operations */
        struct edit_script *next;
} edit_script;

typedef struct edit_script_list {
        int   offset1, offset2;
        int   len1, len2; 
        int    score;
        struct edit_script *script;
        struct edit_script_list *next_script;
} edit_script_list;

struct edit {
        struct edit *link;              /* previous edit instruction */
        char   type[2];
        int    accumulation;
        char   op;                      /* INSERT, DELETE or INTRON */
        int    line1;                   /* line number in file1 */
        int    line2;                   /* line number in file2 */
};

typedef  void  *Pointer;

typedef  struct ValNode {
         Pointer      data;
         struct ValNode *next;
} *ValNodePtr;

typedef int signal_t[5][5];

typedef struct spliced {
   int xs, xe, ys, ye, score;
   char type;
   struct spliced *next;
} splice_t;

typedef struct sim4_stats {
   int internal, icoverage, mult, nmatches;
   double fcoverage, marginals;
} sim4_stats_t;

typedef struct sim4_args {
   int ali_flag, poly_flag, acc_flag, reverse, DRANGE, weight, cutoff;
   int set_K, set_C, set_H;
   int W, K, C, X, B, CDS_from, CDS_to;
   char *S;
} sim4_args_t;

#endif
