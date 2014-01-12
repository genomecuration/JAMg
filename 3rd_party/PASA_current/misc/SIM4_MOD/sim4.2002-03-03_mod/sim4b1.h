#ifndef SIM4B1_H
#define SIM4B1_H
/* $Id: sim4b1.h,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $ */

extern uchar *seq1, *seq2;
extern int M, N, encoding[NACHARS];
extern int file_type;
extern coords last_GT, last_CT, last_AG, last_AC;
extern sim4_args_t rs;

Exon *find_previous(Exon *,Exon *);
void   script_flip_list(edit_script_list **);
void link_to_data_list(Pointer,ValNodePtr *,ValNodePtr *);
void ValNodeFreeData(ValNodePtr);
int good_ratio(int);
void flip_list(Exon **,Exon **);
void free_list(Exon *);
void free_table(void);
void free_align(edit_script_list *);
void seq_toupper(uchar *, int, char *); 
void complement_exons(Exon **,int,int);
void print_exons(Exon *);
void print_pipmaker_exons(Exon *,edit_script_list *,char *,int,int,int,int,uchar *,uchar*,int);
void bld_table(uchar *,int,int,int);
struct edit_script_list *SIM4(uchar*,uchar*,int,int,int,int,int,int,int,int*,int*,int*,Exon**,sim4_stats_t *);
#endif /* SIM4B1_H */
