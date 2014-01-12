#ifndef POLY_H
#define POLY_H

#define MIN_EXON 12

#define T_ONLY    1
#define A_ONLY    2
#define BOTH_AT   3

void get_polyAT(uchar *,int,int *,int *,int);

void remove_poly(struct edit_script_list **,Exon *,uchar *,uchar *,int,int *,int *);

#endif
