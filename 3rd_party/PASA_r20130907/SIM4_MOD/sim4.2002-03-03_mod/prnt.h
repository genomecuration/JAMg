#ifndef SIM_PRNT_H
#define SIM_PRNT_H
/* $Id: prnt.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

typedef unsigned int edit_op_t; /* 32 bits */

void print_align_header(SEQ *seq1, SEQ *seq2, argv_scores_t *ds);
void print_align(int score, uchar *seq1, uchar *seq2, int beg1, int end1, int beg2, int end2, int *S);


#endif
