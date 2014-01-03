//Copyright (c) 2003  by  Mihaela Pertea.
#ifndef _BUILD_SMM_H
#define _BUILD_SMM_H


#define PARENT(x) ( (int) ((x)-1)/4 ) /* macro that returns the parent of x in our tree */


typedef struct Orf
{
  char *label;  /* the sequence label */
  char *seq;    /* the orf sequence   */
} tOrf;


typedef struct sNode
{
  int mut_info_pos;  /* the position that has the maximum mutual info     */
  int mut_info_seq;  /* the base that is to be restricted in mut_info_pos */
  int count[INTERVAL][16]; /* count[x][y] where y is {aa ... tt} = the number
                              of times the first base of y occurs in pos x &
                              the last base of y occurs in pos INTERVAL   */
  float prob[4];     /* prob[x] = the probability that base x is in position INTERVAL */
} tNode;





/* function prototypes */
extern void build_model (int, int, int[]);
extern void calc_prob (int *bp_cnt, int sum, double start_prob[], double end_prob[]);

extern int chi2_test (int, int, int[], int);
extern void compute_power4 (int []);
extern void count_base_pairs (int bp_cnt[][16], char *orf, int v);
extern void count_base_pairs_with_restrictions (char *orf, int level, int);
extern void count_window (int count[][16], int start, char *orf);

extern char  Filter  (char);   // added 11 Dec 98

extern double get_mut_info (int *bp_cnt, int sum);
extern void inc_base_pairs (int bp, int *bp_cnt, char stop);

extern void print_node (char print_string [], int level, int start_node, double mut_info, int, int[]);
extern void process_options (int argc, char **argv);
extern tOrf * read_orfs (FILE* Fptr, int *orf_count);

#endif
