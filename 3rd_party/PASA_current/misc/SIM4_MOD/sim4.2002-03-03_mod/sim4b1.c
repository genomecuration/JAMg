/*
 * sim4.c - A program to align a cDNA sequence with a genomic sequence
 *          for the gene.
 */

#ifndef __lint
/*@unused@*/
static const char rcsid[] = 
"$Id: sim4b1.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "psublast.h"

#include "sim4.h"
#include "sim4b1.h"
#include "Xtend1.h"
#include "align.h"
#include "splice.h"
#include "poly.h"

#define  EXTEND_FW (rs.acc_flag?Xextend_fw:extend_fw)
#define  EXTEND_BW (rs.acc_flag?Xextend_bw:extend_bw)

#define  SLIDE_INTRON(x)   (((x)==TRUE)?sync_slide_intron:slide_intron)

uchar      *seq1, *seq2;
int         M, N, encoding[NACHARS];
coords      last_GT, last_CT, last_AG, last_AC;
int         file_type;

sim4_args_t rs;


static int         numMSPs, K, W, X;
static int         G_score, C_score;
static int        *diag_lev;
static Msp_ptr     msp_list, *msp;
static Exon_ptr    exon_list;

static void   merge(Exon **,Exon **); 
static bool   get_sync_flag(Exon *, Exon *, int);
static void   slide_intron(int w, Exon **,uchar *,uchar *);
static void   sync_slide_intron(int w, Exon **,uchar *,uchar *);
static void   wobble(Exon **,Exon **,const char *,const char *,uchar *seq1);
static Exon  *bmatch(uchar *,uchar *,int,int,int,int);
static Exon  *fmatch(uchar *,uchar *,int,int,int,int);
static void   compact_list(Exon **Lblock, Exon **Rblock);
static int    resolve_overlap(Exon *,Exon *,uchar *);
static int    greedy(uchar *,uchar *,int,int,int,int,Exon **, Exon **);
static int    extend_bw(uchar *,uchar *,int,int,int,int,int *,int *);
static int    extend_fw(uchar *,uchar *,int,int,int,int,int *,int *);
static void   pluri_align(int *,int *,Exon *,struct edit_script_list **);
static void   get_stats(Exon *,sim4_stats_t *); 
static int    get_edist(int,int,int,int,uchar *,uchar *);
static int    get_msp_threshold(int len1, int len2);
static int    find_log_entry(long *log4s, int n, int len, int offset);
static Exon  *new_exon(int,int,int,int,int,int,int,Exon *);
static void   add_word(int,int);
static void   extend_hit(int,int,const uchar *const,const uchar * const,int,int,int);
static void   sort_msps(void);
static void   heapify(int,int);
static int    smaller(int,int);
static void   search(uchar *,uchar *,int,int,int);
static int    link_msps(Msp_ptr *msp,int,int,int);
static int    scale(int n);
static void   msp2exons(Msp_ptr *,int,uchar *,uchar *);
static void   free_msps(Msp_ptr **,int *);
static void   exon_cores(uchar*,uchar*,int,int,int,int,int,int,int,int);
static void   relink(Msp_ptr *,int,int,int,int,int,uchar *,uchar *);
static int    dispatch_find_ends(int,int,int *,int *,edit_script_list *,int,int,int);
static int    find_ends(edit_script_list *,int);
static bool   get_match_quality(Exon *,Exon *,sim4_stats_t *,int);
static int    check_consistency_intron_ori(Exon *,int,char *);

Exon  *find_previous(Exon *,Exon *);
void   script_flip_list(edit_script_list **);


#ifdef DEBUG
static void   debug_print_exons(Exon *, const char *);
#endif

/*  Not currently used: */
#ifdef AUXUTILS 
   static void   remove_polyA_tails(Exon *,uchar *,uchar *,int);
   static void   find_introns(Exon *, Intron **);
   static void   print_introns(Intron *);
#endif

/* seq1 = genomic  DNA (text); seq2 = cDNA */
struct edit_script_list *SIM4(uchar *in_seq1, uchar *in_seq2, 
                         int in_M, int in_N, int in_W, 
                         int in_X, int in_K, int in_C, int in_H,
                         int *dist_ptr, int *pT, int *pA,
                         Exon **Exons, sim4_stats_t *st)
{
  int    cflag, diff, cost, rollbflag, sync_flag;
  int    u, v, I, J;
  bool   good_match;
  Exon   *Lblock, *Rblock=NULL, *tmp_block, *last, *prev, *tmp_block1,
         *tmp_Lblock=NULL, *tmp_Rblock=NULL, *new;

  struct edit_script_list *Script_head=NULL;
  uchar tmp[50];
  coords *sig;

  seq1 = in_seq1;
  seq2 = in_seq2;
  M = in_M;
  N = in_N;
  W = in_W;
  X = in_X;

  if (M<=0 || in_N<=0) {
      *Exons = NULL; return NULL;
  }

  if (rs.acc_flag) {
      last_AG.pos1 = last_AG.pos2 = last_AC.pos1 = last_AC.pos2 = 0;
      last_GT.pos1 = last_GT.pos2 = last_CT.pos1 = last_CT.pos2 = 0;
  }

  /* Compute the distance between two sequences A and B */

  *dist_ptr = 0;

  exon_cores(seq1-1, seq2-1, M, N, 1, 1, 0, W, in_K, PERM); 

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
     if (tmp_block->next_exon==NULL)
         Rblock = tmp_block;
     tmp_block = tmp_block->next_exon;
  }

  if (Lblock && 
      ((Lblock->from1>50000 && Lblock->from2>100) || 
       ((M-Rblock->to1>50000) && (N-Rblock->to2>100)))) {
      free_list(exon_list);
      relink(msp,numMSPs,(in_H>0) ? in_H:DEFAULT_RELINK_H,1,1,0,seq1,seq2);
      tmp_block = Lblock = exon_list;
      while (tmp_block) {
         if (tmp_block->next_exon==NULL)
             Rblock = tmp_block;
         tmp_block = tmp_block->next_exon;
      }
  }
  free_msps(&msp, &numMSPs);

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
     if (tmp_block->next_exon==NULL) 
         Rblock = tmp_block;
     tmp_block = tmp_block->next_exon; 
  }

  /* enclose the current path in the (0,0,0,0) and (M+1,N+1,0,0) brackets */

  Lblock = new_exon (0,0,0,0,0,0,0,Lblock);
  if (Rblock == NULL)
      Rblock = Lblock;
  Rblock->next_exon = new_exon (M+1,N+1,0,0,0,0,0,NULL); 

  /* compute current statistics */
  good_match = get_match_quality(Lblock, Rblock, st, N);

#ifdef DEBUG
  debug_print_exons(Lblock, "LSIS");
#endif

  tmp_block = Lblock;
  while ((tmp_block1 = tmp_block->next_exon)!=NULL) {
     rollbflag = 0;
     diff = (int)(tmp_block1->from2-tmp_block->to2-1);
     if (diff) {
       if (diff<0) {
         int best_u;
         
         best_u = resolve_overlap(tmp_block,tmp_block1,seq1);

         tmp_block1->from1 += best_u+1-tmp_block1->from2;
         tmp_block1->from2 = best_u+1;
         if (((u=tmp_block1->to2-tmp_block1->from2+1)<=0) || (u<8) ||
             ((v=tmp_block1->to1-tmp_block1->from1+1)<=0) || (v<8)) { 
             /* remove exon associated with tmp_block1 */ 
             tmp_block->next_exon = tmp_block1->next_exon;
             tmp_block->flag = tmp_block1->flag; 
             rollbflag = 1;
             free(tmp_block1);
             tmp_block1 = NULL; /* not necessary, just to keep it 'clean'*/
         }
             
         tmp_block->to1 -= tmp_block->to2-best_u; 
         tmp_block->to2 = best_u;
         if (((u=tmp_block->to2-tmp_block->from2+1)<=0) || (u<8) || 
             ((v=tmp_block->to1-tmp_block->from1+1)<=0) || (v<8)) {
             
             /* remove exon defined by tmp_block */
             prev = find_previous(Lblock,tmp_block);
             assert (prev!=NULL);
             prev->next_exon = tmp_block->next_exon;
             prev->flag = tmp_block->flag; 
             if (u>0) rollbflag = 1;
             free(tmp_block);
             tmp_block = prev;
         }

         if (tmp_block->to1)
             tmp_block->length = tmp_block->to2-tmp_block->from2+1;
         if (tmp_block1 && tmp_block1->to1)
             tmp_block1->length = tmp_block1->to2-tmp_block1->from2+1;

       } else {
         /* bridge the gap */
         cflag = (tmp_block1->to2 && tmp_block->to2) ? 0 : 1;
         if (diff && (tmp_block1->from1-tmp_block->to1-1>0)) {
             if (!cflag) {
              if (diff<=MAX_GRINIT) {
                cost = greedy(seq2+tmp_block->to2,
                              seq1+tmp_block->to1, 
                              diff,
                              tmp_block1->from1-tmp_block->to1-1,
                              tmp_block->to2,tmp_block->to1,
                              &tmp_Lblock, &tmp_Rblock);
              } else cost = max(W,(int)(P*diff+1))+1;

              if (cost>max(W,(int)(P*diff+1))) {
                  if (!tmp_block->flag && !tmp_block1->flag) {
                      exon_cores(seq1+tmp_block->to1-1,
                                 seq2+tmp_block->to2-1,
                                 tmp_block1->from1-tmp_block->to1-1,
                                 diff,
                                 tmp_block->to1+1,
                                 tmp_block->to2+1,
                                 1,
                                 min(8,W),
                                 in_C,
/*                               (min(8,W)==W) ? PERM : TEMP); */
                                 TEMP);
                      tmp_Lblock = tmp_Rblock = exon_list;
                      while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL)) 
                         tmp_Rblock = tmp_Rblock->next_exon;

                      if ((!tmp_Lblock && tmp_block1->from1-tmp_block->to1>50000) ||
                           (tmp_Lblock && (tmp_Lblock->from2-tmp_block->to2>100) && 
                           (tmp_Lblock->from1-tmp_block->from1>50000)) || 
                          (tmp_Lblock && (tmp_block1->from2-tmp_Rblock->to2>100) &&
                           (tmp_block1->from1-tmp_Rblock->from1>50000))) {
                         /* possible large intron; increase the score weight */
                         free_list(tmp_Lblock); 
                         relink(msp, numMSPs,
                                (in_H>0) ? in_H:DEFAULT_RELINK_H,
                                tmp_block->to1+1,
                                tmp_block->to2+1,
                                1, seq1, seq2); 
                         tmp_Lblock = tmp_Rblock = exon_list;
                         while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                            tmp_Rblock = tmp_Rblock->next_exon;
                      }    
                      free_msps(&msp, &numMSPs);

                      if (tmp_Lblock) rollbflag = 1;
                      else rollbflag = 0;   /* already 0 */
                  } else 
                        tmp_Lblock = tmp_Rblock = NULL;
              }
             } else if (tmp_block1->to1) {
                   /* start of seq; find last_AG, last_AC */
                   if (rs.acc_flag) {
                     for (v=tmp_block1->from1-1; v<=tmp_block1->to1-3; v++) 
                        if (!strncmp((char *)(seq1+v-2),"AG",(size_t)2)) {
                          last_AG.pos1 = v+1;
                          last_AG.pos2 = tmp_block1->from2+
                                         (v-tmp_block1->from1)+1; 
                          break;
                        } 
                     for (v=tmp_block1->from1-1; v<=tmp_block1->to1-3; v++)
                        if (!strncmp((char *)(seq1+v-2),"AC",(size_t)2)) {
                          last_AC.pos1 = v+1;
                          last_AC.pos2 = tmp_block1->from2+
                                         (v-tmp_block1->from1)+1;
                          break;
                        }
                   } /* end acc_flag */       

                   diff = (int)(min(diff,(int)(MAX_GRINIT/2)));
                   u = min(4*diff,tmp_block1->from1-tmp_block->to1-1); 
                   cost = EXTEND_BW(seq2+tmp_block->to2+
                                    (tmp_block1->from2-tmp_block->to2-1)-diff,
                                    seq1+tmp_block->to1+
                                    (tmp_block1->from1-tmp_block->to1-1)-u,
                                    (int)diff, u,
                                    tmp_block->to2+
                                    (tmp_block1->from2-tmp_block->to2-1)-diff,
                                    tmp_block->to1+
                                    (tmp_block1->from1-tmp_block->to1-1)-u,
                                    &I, &J);
                   if ((good_match==FALSE) || tmp_block->flag || (J==0) || (I==0)) {
                       tmp_block1->from2 = I+1;
                       tmp_block1->from1 = J+1;
                       tmp_block1->edist += cost;
                       tmp_block1->length = tmp_block1->to2-tmp_block1->from2+1;
                   }

                   /* use blast if marginal gap still exists, and this is first scan */
                   if (!(diff=(int)(tmp_block1->from2-tmp_block->to2-1)) ||
                       tmp_block->flag) {
                       /* blast-treated region or no gap */
                       tmp_Rblock = tmp_Lblock = NULL;
                   } else {
                      exon_cores(seq1+tmp_block->to1-1,
                                 seq2+tmp_block->to2-1,
                                 tmp_block1->from1-tmp_block->to1-1,
                                 diff,
                                 tmp_block->to1+1,
                                 tmp_block->to2+1,
                                 1,
                                 min(10,W),
                                 in_C,
/*                               (min(10,W)==W) ? PERM : TEMP); */
                                 TEMP);
                      tmp_block -> flag = 1;
                      tmp_Lblock = tmp_Rblock = exon_list;
                      while (tmp_Rblock && tmp_Rblock->next_exon)
                         tmp_Rblock = tmp_Rblock->next_exon;

                      if ((!tmp_Lblock && tmp_block1->from1-tmp_block->to1>50000) ||
                           (tmp_Lblock && (tmp_Lblock->from2-tmp_block->to2>100) &&
                           (tmp_Lblock->from1-tmp_block->from1>50000)) ||
                          (tmp_Lblock && (tmp_block1->from2-tmp_Rblock->to2>100) &&
                           (tmp_block1->from1-tmp_Rblock->from1>50000))) {
                         /* possible large intron; increase the score weight */
                         free_list(tmp_Lblock);
                         relink(msp, numMSPs,
                                (in_H>0) ? in_H:DEFAULT_RELINK_H,
                                tmp_block->to1+1,
                                tmp_block->to2+1,
                                1,seq1,seq2);
                         tmp_Lblock = tmp_Rblock = exon_list;
                         while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                            tmp_Rblock = tmp_Rblock->next_exon;
                      }
                      free_msps(&msp, &numMSPs);

                      if (tmp_Lblock) rollbflag = 1;
                      else {
                         tmp_block1->from2 = I+1;
                         tmp_block1->from1 = J+1;
                         tmp_block1->edist += cost;
                         tmp_block1->length = tmp_block1->to2-tmp_block1->from2+1;
                      }
                   }

             } else {
                   if (rs.acc_flag) {
                     for (v=tmp_block->to1; v>=tmp_block->from1; v--)
                        if (!strncmp((char *)(seq1+v),"GT",(size_t)2)) {
                          last_GT.pos1 = v;
                          last_GT.pos2 = tmp_block->to2-(tmp_block->to1-v);
                          break;
                        }
                     for (v=tmp_block->to1; v>=tmp_block->from1; v--)
                        if (!strncmp((char *)(seq1+v),"CT",(size_t)2)) {
                          last_CT.pos1 = v;
                          last_CT.pos2 = tmp_block->to2-(tmp_block->to1-v);
                          break;
                        }
                   }

                   diff = (int)(min(diff,(int)(MAX_GRINIT/2)));
                   cost = EXTEND_FW(seq2+tmp_block->to2,
                                    seq1+tmp_block->to1,
                                    diff,
                                    min(4*diff,tmp_block1->from1-tmp_block->to1-1),
                                    tmp_block->to2,tmp_block->to1,
                                    &I, &J);
                   if ((good_match==FALSE) || tmp_block1->flag || (I==M) || (J==N)) {
                       if (tmp_block->to1) {
                          tmp_block->to2 = I;
                          tmp_block->to1 = J;
                          tmp_block->edist += cost;
                          tmp_block->length = tmp_block->to2-tmp_block->from2+1;
                          tmp_Rblock = tmp_Lblock = NULL;
                       } else 
                          /* special case: no initial exon */
                          tmp_Lblock = tmp_Rblock = NULL;
                   }
                   /* use blast if marginal gap still exists, and this is first scan */
                   if (!(diff=(int)(tmp_block1->from2-tmp_block->to2-1)) ||
                        tmp_block1->flag) {
                       /* blast-treated region or no gap */
                       tmp_Rblock = tmp_Lblock = NULL;
                   } else {
                      exon_cores(seq1+tmp_block->to1-1,
                                 seq2+tmp_block->to2-1,
                                 tmp_block1->from1-tmp_block->to1-1,
                                 diff,
                                 tmp_block->to1+1,
                                 tmp_block->to2+1,
                                 1,
                                 min(10,W),
                                 in_C,
/*                               (min(10,W)==W) ? PERM : TEMP); */
                                 TEMP);
                      tmp_Lblock = tmp_Rblock = exon_list;
                      while (tmp_Rblock && tmp_Rblock->next_exon)
                         tmp_Rblock = tmp_Rblock->next_exon;

                      if ((!tmp_Lblock && tmp_block1->from1-tmp_block->to1>50000) ||
                           (tmp_Lblock && (tmp_Lblock->from2-tmp_block->to2>100) &&
                           (tmp_Lblock->from1-tmp_block->from1>50000)) ||
                          (tmp_Lblock && (tmp_block1->from2-tmp_Rblock->to2>100) &&
                           (tmp_block1->from1-tmp_Rblock->from1>50000))) {
                         /* possible large intron; increase the score weight */
                         free_list(tmp_Lblock);
                         relink(msp, numMSPs,
                                (in_H>0) ? in_H:DEFAULT_RELINK_H,
                                tmp_block->to1+1,
                                tmp_block->to2+1,
                                1,seq1,seq2);
                         tmp_Lblock = tmp_Rblock = exon_list;
                         while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                            tmp_Rblock = tmp_Rblock->next_exon;
                      }
                      free_msps(&msp, &numMSPs); 

                      tmp_block1->flag = 1;
                      if (tmp_Lblock) rollbflag = 1;
                      else {
                         if (tmp_block->to1) {
                            tmp_block->to2 = I;
                            tmp_block->to1 = J;
                            tmp_block->edist += cost;
                            tmp_block->length = tmp_block->to2-tmp_block->from2+1;
                            tmp_Rblock = tmp_Lblock = NULL;
                         } else
                            /* special case: no initial exon */
                            tmp_Lblock = tmp_Rblock = NULL;
                      }
                   }
             } 
           } else if (diff) {
                  tmp_Rblock = tmp_Lblock = NULL;
           }
   
           /* merge block in the exon list; make connections 
              to the previous list of blocks; maintain 
              increasing order */

           if (tmp_Lblock) {       
               tmp_block->next_exon = tmp_Lblock;
               tmp_Rblock->next_exon = tmp_block1;
               merge(&tmp_block,&tmp_block1);
           }
       }
    }  /* diff!=0 */
    if (!rollbflag) tmp_block = tmp_block1;
  }

  /* just printing ... */
#ifdef DEBUG
  debug_print_exons(Lblock, "EXTENSIONS");
#endif

  /* compaction step; note: it resets the right end of the list to   */ 
  /* the last item in the block list                                 */

  compact_list(&(Lblock->next_exon), &Rblock);

  /* just printing ... */
#ifdef DEBUG
  debug_print_exons(Lblock, "NORMALIZATION");
#endif

  /* eliminate marginal small blocks at the start of the sequence;   */
  /* resets the empty alignment to one block (Lblock) only           */

  tmp_block = Lblock->next_exon;

  while ((tmp_block!=NULL) && (tmp_block->length<W) && tmp_block->to1) {
          tmp_block1 = tmp_block; /* free */
          tmp_block = tmp_block->next_exon;
          free(tmp_block1); /* free */
  }
  Lblock->next_exon = tmp_block;

  /* eliminate marginal small blocks at the end of the sequence      */

  last = Lblock->next_exon;
  tmp_block = last;
  while (tmp_block!=NULL) {
     if (tmp_block->length>=W)
         last = tmp_block;
     tmp_block = tmp_block->next_exon;
  }
  if (last && last->to1)
      last->next_exon = Rblock->next_exon;
  Rblock = last;

  /* if high accuracy requirement, adjust boundaries of marginal exons */
  if (rs.acc_flag) {
      tmp_block = Lblock->next_exon;

      /* condition for non-signal */
      if (tmp_block && tmp_block->to1 &&
          (strncmp((char *)(seq1+tmp_block->from1-3), END_SIG, (size_t)2) || 
          (tmp_block->from2!=1))) {
          sig = (G_score>=abs(C_score)) ? &last_AG : &last_AC; 
          if (sig->pos1 && (sig->pos2<=20)) {
              /* generated in extend_bw */
	      assert(sig->pos2 > 1);
              (void)strcpy((char *)tmp,END_SIG);
              (void)strncpy((char *)(tmp+2),(char *)seq2,(size_t)sig->pos2-1);
              (void)strcpy((char *)(tmp+sig->pos2+1), START_SIG);
              new = bmatch(seq1,tmp,tmp_block->from1-3,sig->pos2+3,1,1); 
              if (new) {
                  Lblock->next_exon->from1 = sig->pos1;
                  Lblock->next_exon->from2 = sig->pos2;
                  Lblock->next_exon->length -= sig->pos2-1;
                  new->next_exon = Lblock->next_exon;
                  new->ori = (G_score>=abs(C_score)) ? 'G' : 'C';
                  Lblock->next_exon = new; 
              }
          }
      }
      while (tmp_block && tmp_block->next_exon && tmp_block->next_exon->to1) 
             tmp_block = tmp_block->next_exon;
      if (tmp_block && tmp_block->to1 &&
          (strncmp((char *)(seq1+tmp_block->to1),START_SIG,(size_t)2) || (tmp_block->to2!=N))) {   
          sig = (G_score>=abs(C_score)) ? &last_GT : &last_CT;
          if (sig->pos1 && (N-sig->pos2<=20)) {
	      assert(N-sig->pos2 >= 0);
              (void)strcpy((char *)tmp,END_SIG);
              (void)strncpy((char *)(tmp+2),(char *)(seq2+sig->pos2),
			    (size_t)N-sig->pos2);
              (void)strcpy((char *)(tmp+N-sig->pos2+2),START_SIG);
              new = fmatch(seq1+sig->pos1-1,tmp,
                           M-sig->pos1+1,N-sig->pos2+4,
                           sig->pos1-1,sig->pos2+1);
              if (new) {
                  tmp_block->to1 = sig->pos1;
                  tmp_block->to2 = sig->pos2;
                  new->next_exon = tmp_block->next_exon;
                  tmp_block->next_exon = new;
                  tmp_block->ori = (G_score>=abs(C_score)) ? 'G' : 'C';
              }
          }
      } 
  }

  /* Slide exon boundaries for optimal intron signals */
  sync_flag = get_sync_flag(Lblock, Rblock, 6);
  SLIDE_INTRON(sync_flag)(6,&Lblock,seq1,seq2);

  /* decreasingly; script will be in reverse order */
  flip_list(&Lblock, &Rblock); 
  pluri_align(dist_ptr,&(st->nmatches),Lblock,&Script_head); 
  flip_list(&Lblock, &Rblock);      /* increasingly */

  if (rs.poly_flag) 
      remove_poly(&Script_head,Lblock,seq1,seq2,N,pT,pA);
  else
      *pT = *pA = 0;

  get_stats(Lblock, st);

  *Exons = Lblock->next_exon;
  free(Lblock);
  if (!rs.ali_flag) {
     free_align(Script_head);
     return NULL;
  } else
     return Script_head;

}


struct hash_node {
        int ecode;             /* integer encoding of the word */
        int pos;                /* positions where word hits query sequence
 */
        struct hash_node *link; /* next word with same last 7.5 letters */
};

#define HASH_SIZE 32767        /* 2**15 - 1 */
#define GEN_LOG4_ENTRIES 45
#define CDNA_LOG4_ENTRIES 25
static struct hash_node *phashtab[HASH_SIZE+1];
static struct hash_node **hashtab;
static int mask;
static int *next_pos, *pnext_pos;

/* The log4 arrays were computed to mimick the behaviour of the log formula
   for computing the msp threshold in exon_cores(). For genomic_log4s,
   entry i stores the value for the length of a genomic sequence
   for which the contribution to the msp threshold is i/2, i.e.:
               1.4*log_4(3/4*len1) = i/2;  
   Similarly, cDNA_log4s entries store lengths of the cDNA sequence for which
   the contribution to the msp threshold is i/2, i.e.:
               1.4*log_4(len2) = i/2;
   Both arrays are sorted in increasing order, and can be searched with 
   binary search.
*/
static long genomic_log4s[]= {1, 2, 3, 5, 9, 15, 26, 42, 70, 114, \
            188, 309, 507, 832, 1365, 1365, 2240, 2240, 3675, 6029,\
            9892, 16231, 26629, 43690, 71681, \
            117606, 192953, 316573, 519392, 852152,
            1398101, 2293823, 3763409, 6174516, 10130347, \
            16620564, 27268873, 44739242, 73402365, 120429110, \
            197584514, 324171126, 531858072, 872603963, 1431655765
            };
static long cDNA_log4s[]= {1, 1, 2, 4, 7, 11, 19, 32, 52, 86, \
            141, 231, 380, 624, 1024, 1680, 2756, 4522, 7419, 12173, \
            19972, 32768, 53761, 88204, 144715
            };


static int get_msp_threshold(int len1, int len2)
{
    int i, j;

    i = find_log_entry(genomic_log4s, GEN_LOG4_ENTRIES, len1, 0);
    j = find_log_entry(cDNA_log4s, CDNA_LOG4_ENTRIES, len2, 0);

    if (!(i % 2)) return (int)(i/2+j/2);
    else if (!(j % 2)) return (int)(i/2+j/2);
    else return (int)(i/2+j/2+1);
}

static int find_log_entry(long *log4s, int n, int len, int offset)
{
    int a;

    a = n/2;  
    if ((len<log4s[a]) && (!a || (len>=log4s[a-1]))) 
                return max(0,(a-1))+offset;
    else if ((len>=log4s[a]) && ((a==n-1) || (len<log4s[a+1]))) 
                return min(n-1,(a+1))+offset;
    else if (len<log4s[a]) 
                return find_log_entry(log4s,a-1,len, offset);   
    else if (len>log4s[a])
                return find_log_entry(log4s+a+1,n-a-1,len, offset+a+1);
    return -1;
}

/* --------------------   exon_cores()   --------------------- */

static void  exon_cores(uchar *s1, uchar *s2, int len1, int len2, int offset1, int offset2, int flag, int in_W, int in_K, int type)
{
    int       i, W, last_msp, lower, upper;
    int      *allocated;
    Exon     *tmp_block;


    if (in_K<=0) {
        /* compute expected length of longest exact match .. */
        /* K = (int) (log(.75*(double)len1)+log((double)len2))/log(4.0); */
        /* .. and throw in a fudge factor */
        /* K *= 1.4; */

        K = get_msp_threshold(len1, len2);
        if (K>=0) K--; /* compensate for the rounding in the log formula */
/*      commented this to avoid fragmentation 
        if (flag) K = min(K, DEFAULT_C); second pass 
*/
    } else 
        K = in_K;

    numMSPs = 0;
    exon_list = NULL;

    allocated = ckalloc((len1+len2+1)*sizeof(int));
    lower = ((file_type==EST_GEN) || (file_type==GEN_EST && type==TEMP)) 
             ? -len1 : -len2;   
    upper = ((file_type==EST_GEN) || (file_type==GEN_EST && type==TEMP)) 
             ?  len2 :  len1;     
    diag_lev = allocated - lower;
    for (i=lower; i<=upper; ++i) diag_lev[i]=0;

    W = min(in_W,len2);
    switch (file_type) {
       case EST_GEN:  bld_table(s2,len2,W, type);
                      /* use longer sequence for permanent tables */
                      search(s1,s2,len1,len2,W);
                      break;
       case GEN_EST:  if (type!=TEMP) {
                      uchar *aux; int   auxi;

                      aux = s1; s1 = s2; s2 = aux; 
                      auxi = len1; len1 = len2; len2 = auxi;
                      }
                      bld_table(s2,len2,W, type);
                      /* use longer sequence for permanent tables */
                      search(s1,s2,len1,len2,W); 
                      if (type!=TEMP) {
                          register int   auxi; 
                          uchar *aux;
                          Msp_ptr mp; 

                          /* change s1 and s2 back */
                          aux = s1; s1 = s2; s2 = aux;
                          auxi = len1; len1 = len2; len2 = auxi;

                          for (mp=msp_list, i=0; i<numMSPs; i++) {
		      	      auxi = mp->pos1;
			      mp->pos1 = mp->pos2;
			      mp->pos2 = auxi;
			      mp = mp->next_msp;
                          }
                      }
                      break;
       default:       fatal("sim4b1.c: Invalid file type code.");
    }
  
    free(allocated);
    if (type==TEMP) {
        register struct hash_node *hptr, *tptr;
        register int    hval;
        
        free(next_pos);
        for (hval=0; hval<HASH_SIZE+1; hval++) {
             hptr = hashtab[hval];
             while (hptr) {
                tptr = hptr;
                hptr = hptr->link;
                free(tptr);
             }
        }    
        free(hashtab);
    }   

    msp = (Msp_ptr *) ckalloc(numMSPs*sizeof(Msp_ptr));
    { Msp_ptr mp = msp_list;
      for (i = 0; i < numMSPs; ++i) {
                msp[i] = mp;
		mp = mp->next_msp;
      }
    }

    sort_msps();    /* sort in order of mp->pos2, in the shorter seq */

    /* organize Blast hits (MSPs) into exons */
    last_msp = link_msps(msp, numMSPs, DEFAULT_WEIGHT, LINK);

#ifdef DEBUG
    for (i = last_msp; i >= 0; i = msp[i]->prev)
         (void)printf("%d-%d\n", msp[i]->pos1, msp[i]->pos1 + msp[i]->len - 1);
#endif

    msp2exons(msp,last_msp,s1,s2); 

    /* now free msp[]? No - may need to re-link */
    /*
    for (i=0; i<numMSPs; ++i)
             free(msp[i]);
    free(msp);
    */

    tmp_block = exon_list;
    while (tmp_block!=NULL) {
       tmp_block->length = tmp_block->to2-tmp_block->from2+1;
       tmp_block->to1 += offset1;
       tmp_block->from1 += offset1;
       tmp_block->to2 += offset2;
       tmp_block->from2 += offset2;
       tmp_block->flag = flag;

       tmp_block = tmp_block->next_exon;
    }

    return ;
}

static void relink(Msp_ptr *in_msp, int in_numMSPs, int H, int offset1, int offset2, int flag, uchar *s1, uchar *s2)
{
    int last_msp;
    Exon *tmp_block;

    exon_list = NULL;

    last_msp = link_msps(in_msp, in_numMSPs, H, RELINK);

    msp2exons(in_msp,last_msp,s1,s2); 

    tmp_block = exon_list;
    while (tmp_block!=NULL) {
       tmp_block->length = tmp_block->to2-tmp_block->from2+1;
       tmp_block->to1 += offset1;
       tmp_block->from1 += offset1;
       tmp_block->to2 += offset2;
       tmp_block->from2 += offset2;
       tmp_block->flag = flag;

       tmp_block = tmp_block->next_exon;
    }

    return ;
}

static void free_msps(Msp_ptr **in_msp, int *in_numMSPs)
{
     int i;

     for (i=0; i<*in_numMSPs; ++i) free((*in_msp)[i]);
     free(*in_msp); *in_msp = NULL; *in_numMSPs = 0;
}

static int scale(int n)
{
    return (n<=100000) ? n : (100000+(int)(10*log((double)(n-100000))));
}

static int link_msps(Msp_ptr *msp, int numMSPs, int H, int flag)
{
	int i, j, f1, f2, best, diag, diff_diag, best_sc, try;
	

	for (best = -1, best_sc = MININT, i = 0; i < numMSPs; ++i) {
		f1 = msp[i]->pos1;      /* start position in seq1 */
		f2 = msp[i]->pos2;      /* start position in seq2 */
		diag = f1 - f2;
		msp[i]->prev = -1;
		msp[i]->Score = 0;
		for (j = 0; j < i; ++j) {
                     int var_L = 
                     ((msp[i]->pos2+msp[i]->len-msp[j]->pos2-msp[j]->len>2*W) &&
                      (msp[i]->pos2-msp[j]->pos2>2*W)) ? 2*L : L; 
                        
		     diff_diag = diag - msp[j]->pos1 + msp[j]->pos2;
	             if (diff_diag < -rs.DRANGE ||
                         (diff_diag > rs.DRANGE && diff_diag < MIN_INTRON) ||
                         (msp[j]->pos2+msp[j]->len-1-f2>var_L) || 
                         (msp[j]->pos1+msp[j]->len-1-f1>var_L))

				continue;
		     try = msp[j]->Score - ((flag==RELINK) ? 
                           scale(abs(diff_diag)) : abs(diff_diag));

		     if (try > msp[i]->Score) {
		                msp[i]->Score = try;
				msp[i]->prev = j;
		     }
		}
		msp[i]->Score += (H*msp[i]->score);
		if (msp[i]->Score > best_sc) {
			best = i;
			best_sc = msp[i]->Score;
		}
	}

	return best;
}

/* -----------   build table of W-tuples in one of the sequences  ------------*/

void bld_table(uchar *s, int len, int in_W, int type)
{
        int ecode;
        int i, j;
        uchar *t;

        if (type == PERM) {
            mask = (1 << (in_W+in_W-2)) - 1;
            next_pos = pnext_pos;
            hashtab = phashtab;  
            return; 
        }   

        /* perform initializations */
        if (type == INIT) {
 
           /* perform initializations */
           for (i=0; i<NACHARS; ++i) encoding[i]=-1;

           encoding['A'] = 0;
           encoding['C'] = 1;
           encoding['G'] = 2;
           encoding['T'] = 3;
           mask = (1 << (in_W+in_W-2)) - 1;
           /* added +1 in the allocation below */
           next_pos = pnext_pos = (int *)ckalloc((len+1)*sizeof(int));
           hashtab = phashtab;
        } else {
           mask = (1 << (in_W+in_W-2)) - 1;
           next_pos = (int *)ckalloc((len+1)*sizeof(int));
           hashtab = (struct hash_node **)
                      ckalloc((HASH_SIZE+1)*sizeof(struct hash_node *));
           for (i=0; i<=HASH_SIZE; ++i) hashtab[i] = NULL;
        }
    
        /* skip any word containing an N/X  */
        t = s+1;
        for (i=1; (i<=len) && *t; ) {
        restart: 
           ecode = 0L;
           for (j=1; (j<in_W) && (i<=len) && *t; ++j) {
                int tmp = encoding[*t++]; ++i;
                if (tmp<0) goto restart;
                ecode = (ecode << 2) + tmp;
           }
 
           for (; (i<=len) && *t; ) {
                int tmp = encoding[*t++]; i++;
                if (tmp<0) goto restart;   
                ecode = ((ecode & mask) << 2) + tmp;
                add_word(ecode, (int)(t-s-1)); 
           }       
        }

/* previous version: include 'N's in the table  
        t = s;
        ecode = 0L;
        for (i = 1; i < in_W; ++i)
                ecode = (ecode << 2) + encoding[*++t];
        for (; (i<=len) && (*++t); i++) {
                ecode = ((ecode & mask) << 2) + encoding[*t];
                add_word(ecode, (int)(t - s));
        }
 */
}

/* add_word - add a word to the table of critical words */
static void add_word(int ecode, int pos)
{
        struct hash_node *h;
        int hval;

        hval = ecode & HASH_SIZE;
        for (h = hashtab[hval]; h; h = h->link)
                if (h->ecode == ecode)
                        break;
        if (h == NULL) {
                h = (struct hash_node *) ckalloc (sizeof(struct hash_node));
                h->link = hashtab[hval];
                hashtab[hval] = h;
                h->ecode = ecode;
                h->pos = -1;
        }
        next_pos[pos] = h->pos;
        h->pos = pos;
}

/* -----------------------   search the other sequence   ---------------------*/

static void search(uchar *s1, uchar *s2, int len1, int len2, int in_W)
{
        register struct hash_node *h;
        register uchar *t;
        register int ecode, hval;
        int i, j, p;

        t = s1+1;
        for (i=1; (i<=len1) && *t; ) {
        restart:
                ecode = 0L;
                for (j=1; (j<in_W) && (i<=len1) && *t; ++j) {
                    int tmp = encoding[*t++]; ++i;
                    if (tmp<0) goto restart;
                    ecode = (ecode << 2) + tmp;
                }
                for (; (i<=len1) && *t; ) {
                    int tmp = encoding[*t++]; ++i;
                    if (tmp<0) goto restart;
                    ecode = ((ecode & mask) << 2) + tmp;
                    hval = ecode & HASH_SIZE;
                    for (h = hashtab[hval]; h; h = h->link)
                         if (h->ecode == ecode) {
                             for (p = h->pos; p >= 0; p = next_pos[p])
                                  extend_hit((int)(t-s1-1),p,s1,s2,len1,len2, in_W);
                             break;
                         }
                }
        }     

/* previous version: allow 'N's in the table -------- 
        t = s1;
        ecode = 0L;
        for (i = 1; i < in_W; ++i)
                ecode = (ecode << 2) + encoding[*++t];
        for (; (i<=len1) && (*++t); i++) {
                ecode = ((ecode & mask) << 2) + encoding[*t];
                hval = ecode & HASH_SIZE;
                for (h = hashtab[hval]; h; h = h->link)
                        if (h->ecode == ecode) {
                                for (p = h->pos; p >= 0; p = next_pos[p])
                                        extend_hit((int)(t-s1),p,s1,s2,len1,len2, in_W);
                                break;
                        }
        }
---------------------------------------------------*/
}

/* extend_hit - extend a word-sized hit to a longer match */
static void extend_hit(int pos1, int pos2, 
		       const uchar * const s1, const uchar * const s2,
		       int len1, int len2, int in_W)
{
        const uchar *beg2, *beg1, *end1, *q, *s;
        int right_sum, left_sum, sum, diag, score;

        diag = pos2 - pos1;
        if (diag_lev[diag] > pos1)
                return;

        /* extend to the right */
        left_sum = sum = 0;
	q = s1+1+pos1;
        s = s2+1+pos2;
        end1 = q;
        while ((*s != '\0') && (*q != '\0') &&
               (s<=s2+len2) && (q<=s1+len1) && sum >= left_sum - X) {
                sum += ((*s++ == *q++) ? MATCH : MISMATCH);
                if (sum > left_sum) {
                        left_sum = sum;
                        end1 = q;
                }
        }

        /* extend to the left */
        right_sum = sum = 0;
        beg1 = q = (s1+1+pos1) - in_W;
        beg2 = s = (s2+1+pos2) - in_W;
        while ((s>s2+1) && (q>s1+1) && sum >= right_sum - X) {
                sum += ((*(--s) == *(--q)) ? MATCH : MISMATCH);
                if (sum > right_sum) {
                        right_sum = sum;
                        beg2 = s;
                        beg1 = q;
                }
        }

        score = in_W + left_sum + right_sum;
        if (score >= K) {
                Msp_ptr mp = (Msp_ptr)ckalloc(sizeof(*mp));

                mp->len = end1 - beg1;
                mp->score = score;
                mp->pos1 = beg1 - (s1+1);
                mp->pos2 = beg2 - (s2+1);
                mp->next_msp = msp_list;
                msp_list = mp;
                ++numMSPs;
        }

        /*diag_lev[diag] = (end1 - (s1+1)) + in_W; */
        diag_lev[diag] = (end1 - s1) - 1 + in_W;
}

/*  ----------------------------   sort the MSPs  ----------------------------*/

/* sort_msps - order database sequence for printing */
static void sort_msps(void)
{
        int i;
        Msp_ptr mp;

        for (i = (numMSPs/2) - 1; i >= 0; --i)
                heapify(i, (int) numMSPs-1);
        for (i = numMSPs-1; i > 0; --i) {
                mp = msp[0];
                msp[0] = msp[i];
                msp[i] = mp;
                if (i > 1)
                        heapify(0, i-1);
        }
}

/* heapify - core procedure for heapsort */
static void heapify(i, last)
int i, last;
{
        int lim = (last-1)/2, left_son, small_son;
        Msp_ptr mp;

        while (i <= lim) {
                left_son = 2*i + 1;
                if (left_son == last)
                        small_son = left_son;
                else
                        small_son = smaller(left_son, left_son+1);
                if (smaller(i, small_son) == small_son) {
                        mp = msp[i];
                        msp[i] = msp[small_son];
                        msp[small_son] = mp;
                        i = small_son;
                } else
                        break;
        }
}
/* smaller - determine ordering relationship between two MSPs */
static int smaller(i, j)
int i, j;
{
	Msp_ptr ki = msp[i], kj = msp[j];
 
        if (ki->pos2 > kj->pos2) return i;
        if (ki->pos2 < kj->pos2) return j;
        
	return ((ki->pos1 >= kj->pos1) ? i : j);
}

/* ---------------------  organize the MSPs into exons  ---------------------*/

static void msp2exons(Msp_ptr *msp, int last_msp, uchar *s1, uchar *s2)  
{
  Msp_ptr   mp;
  int diag_dist, diff;

  exon_list = NULL; 
  if (last_msp<0) return;

  /* Note: in new_exon, the 'flag' and 'length' fields need not be computed */
  mp = msp[last_msp];
  exon_list = new_exon (mp->pos1, mp->pos2, mp->pos1+mp->len-1, 
                        mp->pos2+mp->len-1, -1, 
                        (mp->len*MATCH-mp->score)/(MATCH-MISMATCH), 0, exon_list);

  last_msp = mp->prev;

  while (last_msp>=0) {
    mp = msp[last_msp]; 
    if (((diag_dist=abs((exon_list->from2-exon_list->from1)-(mp->pos2-mp->pos1)))<=L) 
        && (exon_list->from2-(mp->pos2+mp->len-1))<MAX_INTERNAL_GAP) {
          /* merge with previous exon */
          exon_list->edist += diag_dist;
          exon_list->edist += (mp->len*MATCH-mp->score)/(MATCH-MISMATCH);
          if ((diff=mp->pos2+mp->len-exon_list->from2)>0) {   /* overlap */
             int dist1, dist2;
             dist1 = get_edist(exon_list->from1,mp->pos2+mp->len-diff,
                               exon_list->from1+diff-1,mp->pos2+mp->len-1,s1,s2);
             dist2 = get_edist(mp->pos1+mp->len-diff,mp->pos2+mp->len-diff,
                               mp->pos1+mp->len-1,mp->pos2+mp->len-1,s1,s2);
             exon_list->edist -= max(dist1,dist2);
          } else if (diff<0) {  /* gap */
             exon_list->edist += 0.5*P*(-1)*diff; 
          }
          exon_list->to1 = max(exon_list->to1,mp->pos1+mp->len-1);
          exon_list->to2 = max(exon_list->to2,mp->pos2+mp->len-1);
          exon_list->from1 = min(exon_list->from1,mp->pos1);
          exon_list->from2 = min(exon_list->from2,mp->pos2);
    } else {
          /* new exon */
          exon_list = new_exon (mp->pos1, mp->pos2, mp->pos1+mp->len-1,
                                mp->pos2+mp->len-1, -1,
                                (mp->len*MATCH-mp->score)/(MATCH-MISMATCH),
                                0, exon_list);
    }
    last_msp = mp->prev;
  } 
}

static int get_edist(int f1, int f2, int t1, int t2, uchar *seq1, uchar *seq2)
{
    uchar *s1, *s2, *q1, *q2;
    int dist=0;

    s1 = seq1+f1+1;   /* bc at this stage, the msp pos do not have added +1 */
    s2 = seq2+f2+1;
    q1 = seq1+t1+1;
    q2 = seq2+t2+1;

    while (s1<=q1 && s2<=q2) { dist += (*s1!=*s2); s1++; s2++; } 
    
    return dist;
}

/* ----------------------  print endpoints of exons  --------------------*/

#ifdef AUXUTILS 
static void find_introns(Exon *eleft, Intron **Ilist)
{
  Exon   *tmp_exon, *tmp_exon1;
  Intron *new, *tail;
  int     GTAG_score, CTAC_score;

  *Ilist = tail = NULL;
  if (!eleft) fatal("sim4b1.c: Something wrong in the exon list.\n");

  tmp_exon = eleft->next_exon;
  while ((tmp_exon!=NULL) && (tmp_exon1=tmp_exon->next_exon) &&
         tmp_exon1->to1) {
     
     new = (Intron *)ckalloc(sizeof(Intron)); 
     new->from1 = tmp_exon->to1+1;
     new->to1 = tmp_exon1->from1-1;
     new->from2 = tmp_exon->to2;
     new->to2 = tmp_exon1->from2;
     new->length = new->to1-new->from1+1;
     new->next_intron = NULL;
     if (!tail) *Ilist = new;
     else tail->next_intron = new; 
     tail = new;

     /* find orientation */
     GTAG_score = CTAC_score = 0;
     if (*(seq1+new->from1-1)=='G') GTAG_score++;
     else 
       if (*(seq1+new->from1-1)=='C') CTAC_score++;

     if (*(seq1+new->from1-1)=='T')   
          { GTAG_score++;  CTAC_score++; }

     if (*(seq1+new->to1-1)=='A') 
          { GTAG_score++;  CTAC_score++; }

     if (*(seq1+new->to1-1)=='G')  GTAG_score++;
     else 
       if (*(seq1+new->to1-1)=='C') CTAC_score++;

     if  (GTAG_score>=CTAC_score) 
          new->orientation = '+';
     else new->orientation = 'c';

     tmp_exon = tmp_exon1;
  }
 
}
#endif 


/* should only be called when (file_type==EST_GEN) && (match_ori==BWD) */
void complement_exons(Exon **left, int M, int N)
{ 
  Exon  *tmp_block, *right;
  char   prev, ch;

  prev = 'U'; /* unknown; should trigger error */
  tmp_block = *left;
  while (tmp_block) {
     if (tmp_block->to1) {
         register int aux;
             
         if (tmp_block->next_exon && tmp_block->next_exon->to1) {
               ch = tmp_block->ori;
               tmp_block->ori = prev; 
               switch (ch) {
                 case  'C': prev = 'G'; break;
                 case  'G': prev = 'C'; break;
                 case  'N': prev = 'N'; break;
                 case  'E': prev = 'E'; break;
                 default  : fatal("sim4b1.c: Inconsistency. Check exon orientation at complementation.");
               }
         } else 
               tmp_block->ori = prev;
         aux = tmp_block->from1;
         tmp_block->from1 = M+1-tmp_block->to1;
         tmp_block->to1 = M+1-aux;
         aux = tmp_block->from2;
         tmp_block->from2 = N+1-tmp_block->to2;
         tmp_block->to2 = N+1-aux;
     }
     tmp_block = tmp_block->next_exon;
     if (tmp_block && tmp_block->to1) 
         right = tmp_block;
  }
  flip_list(left,&right);
}

void print_exons(Exon *left)
{
  Exon  *tmp_block, *tmp_block1;

  tmp_block = left;
  while (tmp_block!=NULL) {
     if (tmp_block->to1) {
         
         if (file_type==EST_GEN) 
            (void)fprintf(stdout,"%d-%d  (%d-%d)   %d%%",
                          tmp_block->from2, tmp_block->to2,
                          tmp_block->from1, tmp_block->to1,
                          tmp_block->match);
         else  
            /* file_type==GEN_EST */
            (void)fprintf(stdout,"%d-%d  (%d-%d)   %d%%", 
                           tmp_block->from1, tmp_block->to1,
                           tmp_block->from2, tmp_block->to2,
                           tmp_block->match);
         if (((tmp_block1=tmp_block->next_exon)!=NULL) && tmp_block1->to1)
              switch (tmp_block->ori) {
                 case 'C':  (void)fprintf(stdout," <-\n");
                            break;
                 case 'E':  (void)fprintf(stdout," ==\n");
                            break;
                 case 'G':  (void)fprintf(stdout," ->\n");
                            break;
                 case 'N':  (void)fprintf(stdout," --\n");
                            break;
                 default :  fatal("sim4b1.c: Inconsistency. Check exon orientations.");
              }  
     }        
     tmp_block = tmp_block->next_exon;
  }  
}    

/* to and from are in the original cDNA sequence */
void print_pipmaker_exons(Exon *exons, edit_script_list *aligns, char *gene, 
                          int from, int to, int M, int N, uchar *seq1, uchar *seq2, 
                          int match_ori)
{
   Exon  *tmp_block, *left, *right;
   int From, To, cov, ori;
 
   /* print the first line in the record */
   if ((exons==NULL) || 
       (!exons->to1 && ((exons->next_exon==NULL) || !exons->next_exon))) 
       return;
 
   left = right = tmp_block = (exons->to1) ? exons : exons->next_exon;
   while (tmp_block) {
      if (tmp_block->to1) right = tmp_block; 
      tmp_block = tmp_block->next_exon;
   }
   /* report any inconsistencies between the intron orientations, as well
   as between the mRNA (CDS) match strand and the introns orientations */
   ori = check_consistency_intron_ori(exons, match_ori, gene);

   /* determine the matching coordinates for the CDS */
   if ((from>0) && (to>0) && aligns) {
       cov = dispatch_find_ends(from, to, &From, &To, aligns, M, N, match_ori);
       switch (cov) {
          case OK: if ((match_ori==FWD) && strncmp((char *)(seq1+From-1),"ATG",3))
                      fprintf(stderr, "Warning: No start codon at location %d in the genomic sequence (%s).\n", From, gene);  
                   else if ((match_ori==BWD) && strncmp((char *)(seq1+To-3),"CAT",3)) 
                      fprintf(stderr, "Warning: No (complement) start codon at location %d in the genomic sequence (%s).\n", To-2, gene);
                   if ((match_ori==FWD) && strncmp((char *)(seq1+To-3),"TAA",3) && 
                       strncmp((char *)(seq1+To-3),"TAG",3) && 
                       strncmp((char *)(seq1+To-3),"TGA",3))
                       fprintf(stderr, "Warning: No stop codon at location %d in the genomic sequence (%s).\n", To-2, gene);
                   else if ((match_ori==BWD) && strncmp((char *)(seq1+From-1),"TTA",3) &&
                       strncmp((char *)(seq1+From-1),"CTA",3) && 
                       strncmp((char *)(seq1+From-1),"TCA",3)) 
                       fprintf(stderr, "Warning: No (complement) stop codon at location %d in the genomic sequence (%s).\n", From, gene);
                   break;
     
          case FREE_START: 
                   fprintf(stderr, "Warning: Start of CDS does not match (%s).\n", gene);
                   if ((match_ori==FWD) && strncmp((char *)(seq1+To-3),"TAA",3) && 
                       strncmp((char *)(seq1+To-3),"TAG",3) && 
                       strncmp((char *)(seq1+To-3),"TGA",3))
                       fprintf(stderr, "Warning: No stop codon at location %d in the genomic sequence (%s).\n", To-2, gene);
                   else if ((match_ori==BWD) && strncmp((char *)(seq1+From-1),"TTA",3) &&   
                       strncmp((char *)(seq1+From-1),"CTA",3) &&
                       strncmp((char *)(seq1+From-1),"TCA",3))
                       fprintf(stderr, "Warning: No (complement) stop codon at location %d in the genomic sequence (%s).\n", From, gene);  

                   break;

          case FREE_END:
                   fprintf(stderr, "Warning: End of CDS does not match (%s).\n", gene);
                   if ((match_ori==FWD) && strncmp((char *)(seq1+From-1),"ATG",3))
                      fprintf(stderr, "Warning: No start codon at location %d in the genomic sequence (%s).\n", From, gene);
                   else if ((match_ori==BWD) && strncmp((char *)(seq1+To-3),"CAT",3))
                      fprintf(stderr, "Warning: No (complement) start codon at location %d in the genomic sequence (%s).\n", To-2, gene);

                   break;
    
          case FREE_BOTH_ENDS:
                   fprintf(stderr, "Warning: Start of CDS does not match (%s).\n", gene);
                   fprintf(stderr, "Warning: End of CDS does not match (%s).\n", gene);
                   break;
               
          default: fatal("Unrecognized warning code."); 
       }            
   }

   /* report codon inconsistencies in the cDNA */
   if (to>0 && from>0) {
      if (strncmp((char *)(seq2+from-1),"ATG",3))
         fprintf(stderr, "Warning: No start codon at location %d in the mRNA (%s).\n", from, gene); 
      if (strncmp((char *)(seq2+to-3),"TAA",3) && 
          strncmp((char *)(seq2+to-3),"TAG",3) && 
          strncmp((char *)(seq2+to-3),"TGA",3)) 
         fprintf(stderr, "Warning: No end codon at location %d in the mRNA (%s).\n", to-2, gene);
   }

   printf("%c %d %d %s%s\n",
          (ori==FWD) ? '>':'<', left->from1, right->to1, gene ? gene:"",
          (match_ori==FWD) ? "":" (complement)");
   if ((from>0) && (to>0) && aligns)
          printf("+ %d %d\n", From, To);

   /* now print the exons */
   /* if (match_ori==BWD) flip_list(&left,&right); not accepted by PipMaker */ 
   tmp_block = left;
   while (tmp_block!=NULL) {
      if (tmp_block->to1)
          (void)fprintf(stdout,"%d %d\n",
                           tmp_block->from1, tmp_block->to1);
      tmp_block = tmp_block->next_exon;
   }
   if (match_ori==BWD) flip_list(&left, &right);

   return;
}

static int check_consistency_intron_ori(Exon *exons, int match_ori, char *gene)
{
   Exon *t=exons;
   int numG, numC, numE, numN;
 
   numG = numC = numE = numN = 0;

   if (!t->to1) t = t->next_exon;
   while (t && t->to1) {
     if (t->next_exon  && t->next_exon->to1) {
       switch (t->ori) {
         case 'G': numG++; break;
         case 'C': numC++; break;
         case 'N': numN++; break;
         case 'E': numE++; break;
         default : fatal("sim4b1.c: Unrecognized intron orientation.");
       } 
     }
     t = t->next_exon;
   }
   if (numG && numC) 
    fprintf(stderr, "Warning: Introns reported on both strands (%s).\n", gene);
/*
   Note: a match can be reverse complemented either b/c the complement
   genomic sequence was used as input, while the mRNA was actually transcribed 
   in forward orientation (CT_AC), or b/c the forward strand was given for the
   genomic sequence, but the mRNA was transcribed in the reverse direction 
   (GT_AG); hence there is no relevant test for this  
   else if ((numG && (match_ori==BWD)) || (numC && (match_ori==FWD))) 
    fprintf(stderr, "Warning: Introns orientations inconsistent with the reported match strand (%s).\n", gene);
*/
      
   if (numN) 
    fprintf(stderr, "Warning: Ambiguous intron orientation (%s).\n", gene);
   if (numE) 
    fprintf(stderr, "Warning: Internal gap in the mRNA (%s).\n", gene); 

   return (numG>=numC) ? FWD:BWD;
}


/* from and to are given in the original cDNA sequence */
static int dispatch_find_ends(int from, int to, int *From, int *To, edit_script_list *aligns, int M, int N, int match_ori)
{
   int f1, f2, t1, t2, ot1, ot2, xto, xfrom;
   int free_start, free_end;

   free_start = free_end = 1;

   if (aligns->next_script && (aligns->offset2 > aligns->next_script->offset2))
       script_flip_list(&aligns);
   
   if (match_ori==FWD) {
       xto = to; xfrom = from; 
   } else if (file_type == EST_GEN) {
       xto = to; xfrom = from;
   } else {
       xto = N-from+1; xfrom = N-to+1;
   }

   *From = *To = 0; t1 = t2 = 0;
   while (aligns) {
      ot2 = t2; ot1 = t1;
      f1 = aligns->offset1; f2 = aligns->offset2; 
      t1 = f1+aligns->len1-1; t2 = f2+aligns->len2-1; 
      if (ot2 < xfrom && xfrom < f2) {
         *From = f1; break;
      }
      if (f2 <= xfrom && xfrom <= t2) {
         *From = find_ends(aligns, xfrom); free_start = 0; break;
      }
      aligns = aligns->next_script;
   }
   if (*From == 0) return FREE_BOTH_ENDS;
   
   if (ot2 < xto && xto < f2) *To = ot1; 
   else if (xto <= t2) { *To = find_ends(aligns, xto); free_end = 0; }
   else {
     *To = 0;
     while (aligns && ((aligns=aligns->next_script)!=NULL)) {
        ot2 = t2; ot1 = t1;
        f1 = aligns->offset1; f2 = aligns->offset2; 
        t1 = f1+aligns->len1-1; t2 = f2+aligns->len2-1;
        if (ot2 < xto && xto < f2) { *To = ot1; break; }
        if (t2 < xto) *To = t1;
        else if (f2 <= xto && xto <=t2) {
          *To = find_ends(aligns, xto); free_end = 0; break; 
        }
     }    
     if (*To==0) *To = t1;
   }

   if (*To == 0) { *From = 0; free_start = 1; }
   
   if (*To && *From && match_ori==BWD && file_type==EST_GEN) {
       int aux = M-(*From)+1; *From = M-(*To)+1; *To = aux;
   }

   if (free_start && free_end) return FREE_BOTH_ENDS;
   else if (free_start) return FREE_START;
   else if (free_end) return FREE_END;

   return OK;
}

static int find_ends(edit_script_list *head, int j0)
{
   int i, j, e1, e2;
   edit_script *tp;

   i = head->offset1; e1 = i+head->len1-1;
   j = head->offset2; e2 = j+head->len2-1;

   tp = head->script;
   i--; j--;
   while (i<=e1 && j<=e2 && tp) {
      if (j==j0) return i;
      switch (tp->op_type) {
         case DELETE:     i += tp->num; break;
         case INSERT:     j += tp->num; if (j>=j0) return i; break;
         case SUBSTITUTE: i += tp->num; j += tp->num; 
                          if (j>=j0) return (i-(j-j0)); break;
         default:         fatal("Illegal opcode in script."); 
      }
      tp = tp->next;
   }

   /* not found: failure */
   fatal("Inconsistency in script.");
}
 

#ifdef AUXUTILS
static void print_introns(Intron_ptr intron_list) 
{
  Intron_ptr ep=intron_list;

  (void)printf("\nIntrons {\n\n");
  while (ep!=NULL) {
     (void)printf("genome: %6ld - %-6ld    cDNA: %5ld - %-5ld,  l: %-5d,  o: %c \n",
                   ep->from1, ep->to1, ep->from2, ep->to2, 
                   ep->length,
                   ep->orientation);
     ep = ep->next_intron;
  }
  (void)printf("}\n\n");
}
#endif


static void pluri_align(int *dist_ptr,int *num_matches,Exon *lblock,struct edit_script_list **Aligns)
{
   int    tmpi, di_count, i, end1, end2, diff, ali_dist, nmatches;
   uchar *a, *b;
   Exon  *tmp_block=lblock, *tmp_block1;

   struct edit_script_list *enew;
   struct edit_script *head, *tmp_script, *new, *left, *right, *prev;

   nmatches = 0;

   head = NULL;
   *Aligns = NULL;
   *dist_ptr = ali_dist = 0;

   end1 = M; 
   end2 = N;

   while (((tmp_block1=tmp_block->next_exon)!=NULL) && tmp_block1->to1) {

         if ((diff=tmp_block->from2-tmp_block1->to2-1)!=0) {
           if (tmp_block->to1) {
              enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
              enew->next_script = *Aligns;
              *Aligns = enew;
              (*Aligns)->script = head;
              (*Aligns)->offset1 = tmp_block->from1;
              (*Aligns)->offset2 = tmp_block->from2;
              (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
              (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
              (*Aligns)->score = ali_dist;
              ali_dist = 0;
              head = NULL;
            }
            end1 = tmp_block1->to1;
            end2 = tmp_block1->to2;  

         } else if (((diff=tmp_block->from1-tmp_block1->to1-1)!=0) &&
                      tmp_block->to1) {
              new = (edit_script *) ckalloc(sizeof(edit_script));
              new->op_type = DELETE;
              new->num = diff;
              new->next = head;
              head = new;
         } else if (diff) 
              end1 = tmp_block1->to1;

         diff = align_get_dist(tmp_block1->from1-1, tmp_block1->from2-1,
                               tmp_block1->to1, tmp_block1->to2, 
                               max(1000,.2*(tmp_block1->to2-tmp_block1->from2+1)));

         if (diff<0) {
             (void)printf("The two sequences are not really similar.\n");
             (void)printf("Please try an exact method.\n");
             exit(1);        
         }

#ifdef STATS
         if (diff>P*(tmp_block1->to2-tmp_block1->from2+1))
              (void)printf("Warning: Distance threshold on segment exceeded.\n");
#endif

         align_path(tmp_block1->from1-1, tmp_block1->from2-1, 
                    tmp_block1->to1, tmp_block1->to2, diff, &left, &right);

         Condense_both_Ends(&left, &right, &prev);
  
         if (!tmp_block->to1 && right->op_type == DELETE) {
            /* remove gaps at end of alignment */
            diff -= 0+right->num;         /* subtract GAP_OPEN = 0 */
            tmp_block1->to1 -= right->num;
            end1 -= right->num; 
            if (head && (head->op_type == DELETE)) 
                head->num += right->num;
            free(right); prev->next = NULL; 
            right = prev;
         } 
         if ((!tmp_block1->next_exon || !tmp_block1->next_exon->to1) &&
              left && (left->op_type == DELETE)) {
            diff -= 0+left->num;          /* subtract GAP_OPEN = 0 */
            tmp_block1->from1 += left->num;
            tmp_script = left->next; 
            if (right == left) right = tmp_script;
            free(left); left = tmp_script; 
         }
   
         *dist_ptr += diff;
         ali_dist += diff;

         a = seq1+tmp_block1->from1-1;
         b = seq2+tmp_block1->from2-1;
         tmpi = di_count = 0;
         tmp_script = left;  
         while (tmp_script) {
            switch (tmp_script->op_type) {
                case  DELETE:
                      di_count += tmp_script->num;
                      tmpi += tmp_script->num;
                      a += tmp_script->num;
                      break;
                case  INSERT:
                      di_count += tmp_script->num;
                      tmpi += tmp_script->num;
                      b += tmp_script->num;
                      break;
                case  SUBSTITUTE:
                      for (i=0; i<tmp_script->num; ++i, ++a, ++b)
                        if (*a!=*b) tmpi++; else nmatches++;
                      break;
            }         
            tmp_script = tmp_script->next;
         }
         tmp_block1->alen = (int)((tmp_block1->to1-tmp_block1->from1+1+
                             tmp_block1->to2-tmp_block1->from2+1+di_count)/(double)2);
         tmp_block1->nmatches = tmp_block1->alen - tmpi;
         tmp_block1->match = (int)floor(100*(1-tmpi/(double)tmp_block1->alen));

         right->next = head;
         head = left;
         tmp_block = tmp_block1;
  }


   /* at the beginning of the sequences */
  if (tmp_block1!=NULL) {

    if ((diff=tmp_block->from2-tmp_block1->to2-1)!=0 && (diff!=N)) {
         enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
         enew->next_script = *Aligns;
         *Aligns = enew;
         (*Aligns)->offset1 = tmp_block->from1;
         (*Aligns)->offset2 = tmp_block->from2;
         (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
         (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
         (*Aligns)->script = head;
         (*Aligns)->score = ali_dist;
 
    } else if (diff!=N) {

         /* modified to cut introns at the beginning of the sequence */
         enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
         enew->next_script = *Aligns;
         *Aligns = enew;
         (*Aligns)->offset1 = tmp_block->from1;
         (*Aligns)->offset2 = 1;
         (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
         (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
         (*Aligns)->script = head;
         (*Aligns)->score = ali_dist;
    }
  }
  *num_matches = nmatches;
}

static Exon *new_exon(int f1, int f2, int t1, int t2, int len, int edist, int flag, Exon *next)
{
  Exon *new = (Exon *)ckalloc(sizeof(Exon));

  new->from1 = f1;
  new->from2 = f2;
  new->to1 = t1;
  new->to2 = t2;
  new->length = (len < 0) ? (t2-f2+1) : len;
  new->edist = edist;
  new->flag = flag;
  new->next_exon = next;

  return new;
}


static void get_stats(Exon *lblock, sim4_stats_t *st)
{
   Exon *t, *t1;

#ifdef _STATS
   t = lblock;
   if (!t->next_exon) { /* no alignment found */
       st->marginals = 2.0;
   } else {
     while (t) {
       if ((t1 = t->next_exon)!=NULL) {
        if (!t->to1 && t1->to1)
            st->marginals += (float)(t1->from2-1)/t1->to2;
        else if (t->to1 && !t1->to1)
            st->marginals += (float)(N-t->to2)/
                                (N-t->from2+1);
        else if (!t->to1 && !t1->to1)
            st->marginals = 2.0;
       }
       t = t1;
     }
   }
   st->marginals = st->marginals/2;
#endif

   st->icoverage = 0;  st->internal = 1; st->mult = 0;
   t = lblock->next_exon;
   while (t) {
      st->icoverage += t->length;
      if (t->length) st->mult++;
      t = t->next_exon;
   }
   st->fcoverage = ((float)st->icoverage)/N;

   t = lblock;
   if ((t->next_exon==NULL) || !t->next_exon->to1)                            
      st->internal = 0;
   while (t) {
      if ((t->to1) && ((t1=t->next_exon)!=NULL) && 
          (t1->from2-t->to2-1>0) && t1->to1)
                    st->internal = 0;
      t = t->next_exon;
   }
}

static int resolve_overlap(Exon *tmp_block, Exon *tmp_block1, uchar *seq1)
{          
     int   diff, best_u, l0, l1, u, cost;
     int    GTAG_score, CTAC_score;
     uchar *s1, *s2, *e1;
       
     diff = tmp_block1->from2-tmp_block->to2-1;
     if (diff>=0) return (tmp_block1->from2-1);

     /* resolve overlap using the GT-AG criterion */
     /* u-1 = actual position in the sequence */
       
     l0 = tmp_block->length-diff;
     l1 = tmp_block1->length; 
           
     best_u = u = tmp_block1->from2-1;
     s1 = seq1+tmp_block->to1-(tmp_block->to2-u);
     s2 = seq1-2+tmp_block1->from1+u-tmp_block1->from2;
       
     cost = 0;
     e1 = seq1+tmp_block->to1;
     while (s1<=e1) {
       GTAG_score = CTAC_score = 0;
       GTAG_score += ((char)(*s1)=='G') ? 1 : 0;
       GTAG_score += ((char)(*(s1+1))=='T') ? 1 : 0;
       GTAG_score += ((char)(*s2)=='A') ? 1 : 0;
       GTAG_score += ((char)(*(s2+1))=='G') ? 1 : 0;
     
       if (GTAG_score > abs(cost) && ((l0>=8) || (l1>=8))) {
           cost = GTAG_score;
           best_u = u;
           if (cost == 4) break; 
       }
     
       CTAC_score += ((char)(*s1)=='C') ? 1 : 0;
       CTAC_score += ((char)(*(s1+1))=='T') ? 1 : 0;
       CTAC_score += ((char)(*s2)=='A') ? 1 : 0;
       CTAC_score += ((char)(*(s2+1))=='C') ? 1 : 0;
     
       if (CTAC_score > abs(cost)) {
           cost = -CTAC_score;
           best_u = u;
           if (cost == 4) break;
       }
       
       u++; s1++; s2++;
       l0++; l1--;
     }     
        
     return best_u;
}      


static int  greedy(uchar *s1, uchar *s2, int m, int n, int offset1, int offset2, Exon **lblock, Exon **rblock)
{
  int     col,                    /* column number */
          d,                      /* current distance */
          k,                      /* current diagonal */
          max_d,                  /* bound on size of edit script */
          Cost,
          blower,flower,          /* boundaries for searching diagonals */
          bupper,fupper,
          row,                    /* row number */
          DELTA,                  /* n-m  */
          MAX_D,
          B_ORIGIN, F_ORIGIN;
  int     back, forth;            /* backward and forward limits at exit */

  int     *blast_d, *flast_d,     /* rows containing the last d (at crt step, d-1) */
          *btemp_d, *ftemp_d;     /* rows containing tmp values for the last d */
  int     *min_row, *min_diag,    /* min (b)/ max (f) row (and diagonal) */
          *max_row, *max_diag;    /* reached for cost d=0, ... m.  */

          
  DELTA = n-m;                    
/*max_d = MAX_D = m+1; */         
  max_d = MAX_D = max(W,(int)(P*m+1));
          
  if (DELTA<0) {
      if (m<=min(W,(1+P)*n)) {
          *lblock = *rblock = new_exon(offset2+1,offset1+1,offset2+n,offset1+m,
                                       m,n-m+(int)(P*m+1),0,NULL);
          return m-n+(int)(P*n+1);
      } else {
          *lblock = *rblock = NULL;
          return max(W,(int)(P*m+1))+1;
      }
  }       
  
  F_ORIGIN = MAX_D;
  B_ORIGIN = MAX_D-DELTA;
  for (row=m, col=n; row>0 && col>0 && (s1[row-1]==s2[col-1]); row--,col--)
        /*LINTED empty loop body*/;
  
  if (row == 0) {
      /* hit last row; stop search */
      *lblock = *rblock = new_exon(offset2-m+n+1,offset1+1,offset2+n,
                                   offset1+m,m,0,0,NULL);
      return 0;
  }
          
  
  blast_d = (int *)ckalloc((MAX_D+n+1)*sizeof(int));
  btemp_d = (int *)ckalloc((MAX_D+n+1)*sizeof(int));
          
  for (k=0; k<=MAX_D+n; ++k) { blast_d[k]=m+1; btemp_d[k]=m+1; }
  blast_d[B_ORIGIN+DELTA] = row;
                                   
  blower = B_ORIGIN + DELTA - 1;
  bupper = B_ORIGIN + DELTA + 1;
  
  
  for (row=0; row<n && row<m && (s1[row]==s2[row]); row++)
        /*LINTED empty loop body*/;
  
  if (row == m) {
       /* hit last row; stop search */
       *lblock = *rblock = new_exon(offset2+1,offset1+1,offset2+m,
                                    offset1+m,m,0,0,NULL);
       free(blast_d); free(btemp_d);
   
       return 0;
  }
  
  flast_d = (int *)ckalloc((MAX_D+n+1)*sizeof(int));
  ftemp_d = (int *)ckalloc((MAX_D+n+1)*sizeof(int));
  
  for (k=0; k<=MAX_D+n; ++k) { flast_d[k]=-1; ftemp_d[k]=-1; }
  flast_d[F_ORIGIN] = row;           
  
  flower = F_ORIGIN - 1;
  fupper = F_ORIGIN + 1;
  
  max_row = (int *)ckalloc((MAX_D+1)*sizeof(int));
  min_row = (int *)ckalloc((MAX_D+1)*sizeof(int));
  max_diag = (int *)ckalloc((MAX_D+1)*sizeof(int));
  min_diag = (int *)ckalloc((MAX_D+1)*sizeof(int));
       
  for (d=1; d<=MAX_D; d++) {
       min_row[d] = m+1; max_row[d] = -1;
  }    
  min_row[0] = blast_d[B_ORIGIN+DELTA];
  min_diag[0] = B_ORIGIN+DELTA;
  max_row[0] = flast_d[F_ORIGIN];
  max_diag[0] = F_ORIGIN;
  
  back = forth = -1;
  
  d = 1;
  while (d <= max_d) {             
  
        /* for each relevant diagonal ... */
        for (k = blower; k <= bupper; k++) {
             /* process the next edit instruction */
  
             /* find a d on diagonal k */
             if (k==-d+DELTA+B_ORIGIN) {
  
                      /* move left from the last d-1 on diagonal k+1 */
                      row = blast_d[k+1];   /* INSERT */
             } else if (k==d+DELTA+B_ORIGIN) {
  
                      /* move up from the last d-1 on diagonal k-1 */
                      row = blast_d[k-1]-1;  /* DELETE */
             } else if ((blast_d[k]<=blast_d[k+1]) &&
                        (blast_d[k]-1<=blast_d[k-1])) {
  
                      /* substitution */
                      row = blast_d[k]-1;  /* SUBSTITUTE */
             
             } else if ((blast_d[k-1]<=blast_d[k+1]-1) &&
                        (blast_d[k-1]<=blast_d[k]-1)) {
                      /* move right from the last d-1 on diagonal k-1 */
                      row = blast_d[k-1]-1;  /* DELETE */
             } else  {
                      /* move left from the last d-1 on diagonal k+1 */
                      row = blast_d[k+1];    /* INSERT */
             }
             /* code common to the three cases */
             col = row + k - B_ORIGIN;
                      
             /* slide up the diagonal */
             while (row > 0 && col > 0 && (s1[row-1]==s2[col-1])) {
                             --row;
                             --col;
             }        
             btemp_d[k] = row;
              
/*           if (row == 0 || col == 0) max_d = d;   */
 
        }     /* for k */

        min_row[d] = btemp_d[DELTA+B_ORIGIN];
        min_diag[d] = DELTA+B_ORIGIN;
        for (k=blower; k<=bupper; ++k) {
             blast_d[k] = btemp_d[k]; btemp_d[k] = m+1;
             if (blast_d[k]<min_row[d]) {
                       min_row[d] = blast_d[k];
                       min_diag[d] = k;
             }
        }
 
        /* record cell, if paths overlap with minimum combined cost */
        /* obs: it suffices to search up to Cost=min(d-1,(max_d-d)) */
        for (Cost=0; Cost<d; Cost++) {
             if ((min_row[d]<=max_row[Cost]) &&
                 ((max_d > d+Cost) || (max_d==d+Cost && (forth<0)))) { 
                           max_d = d+Cost;
                           back = d;
                           forth = Cost;
                           break;
             }
        }
             
        --blower; ++bupper;
        
        /* for each relevant diagonal ... */
        for (k = flower; k <= fupper; k++) {
               /* process the next edit instruction */
             
               /* find a d on diagonal k */       
               if (k==-d+F_ORIGIN) {
                        /* move down from the last d-1 on diagonal k+1 */
                        row = flast_d[k+1]+1; /* DELETE */
             
               } else if (k==d+F_ORIGIN) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = flast_d[k-1]; /* INSERT */
                           
               } else if ((flast_d[k]>=flast_d[k+1]) &&
                          (flast_d[k]+1>=flast_d[k-1])) {
        
                        /* substitution */
                        row = flast_d[k]+1;   /* SUBSTITUTE */
               
               } else if ((flast_d[k+1]+1>=flast_d[k-1]) &&
                          (flast_d[k+1]>=flast_d[k])) {
             
                        /* move left from the last d-1 on diagonal k+1 */
                        row = flast_d[k+1]+1;  /* DELETE */
               } else { 
                        /* move right from the last d-1 on diagonal k-1 */
                        row = flast_d[k-1];    /* INSERT */
               }        
               /* code common to the three cases */
               col = row + k - F_ORIGIN;

               /* slide down the diagonal */
               if (row>=0)
               while (row < m && col < n && (s1[row]==s2[col])) {
                                ++row;
                                ++col;
               }        
               ftemp_d[k] = row;
               
/*             if (row == m || col == n) max_d = d;  */
          }     /* for k */    
                   
          max_row[d] = ftemp_d[F_ORIGIN];
          max_diag[d] = F_ORIGIN;
          for (k=flower; k<=fupper; ++k) {
               flast_d[k] = ftemp_d[k]; ftemp_d[k] = -1;
               if (flast_d[k]>max_row[d]) {
                   max_row[d] = flast_d[k];
                   max_diag[d] = k;
               }
          }
                       
          /* record backward and forward limits, if minimum combined
           * cost in overlapping. Note: it suffices to search up to
           * Cost=min(d,(max_d-d)).
           */          
          for (Cost=0; Cost<=d; Cost++) {
             if ((min_row[Cost]<=max_row[d]) &&
                 ((max_d>d+Cost) || (max_d==d+Cost && (forth<0)))) {
                      max_d = d+Cost;
                      back = Cost;
                      forth = d;
                      break;
             }
          }
          --flower; ++fupper;
               
          ++d;  /* for d */
  }
                   
  if (d>MAX_D) {
      *lblock = *rblock = NULL;
                       
      free(blast_d); free(btemp_d);
      free(flast_d); free(ftemp_d);
      free(min_row); free(min_diag);
      free(max_row); free(max_diag);
             
      return d;  
  }                   
/*fin:*/              
  if (m-min_row[back]>=max_row[forth]) {
      *rblock = new_exon(offset2+1+min_row[back]+min_diag[back]-B_ORIGIN,
                         offset1+1+min_row[back],
                         offset2+n,offset1+m,
                         m-min_row[back],back,0,NULL);
      *lblock = new_exon(offset2+1,offset1+1,
                         offset2+min_row[back]+max_diag[forth]-F_ORIGIN,
                         offset1+min_row[back],
                         min_row[back],forth,0,*rblock);
  } else {
      *rblock = new_exon(offset2+1+max_row[forth]+min_diag[back]-B_ORIGIN,
                         offset1+1+max_row[forth],
                         offset2+n,offset1+m,m-max_row[forth],back,0,NULL);
      *lblock = new_exon(offset2+1,offset1+1,
                         offset2+max_row[forth]+max_diag[forth]-F_ORIGIN,
                         offset1+max_row[forth],max_row[forth],forth,0,*rblock);
  }
                      
  free(blast_d); free(btemp_d);
  free(flast_d); free(ftemp_d);
  free(min_row); free(min_diag);
  free(max_row); free(max_diag);

  return back+forth;
}


void flip_list(Exon **left, Exon **right)
{
   Exon   *ep, *ahead, *behind;

   *right = *left;
   ahead = *left;
   ep = NULL;
   while (ahead!=NULL) {
          behind = ep;
          ep = ahead;
          ahead = ahead->next_exon;
          ep->next_exon = behind;
  }
  *left = ep;
}

/* reverse a list of edit script chains */
void script_flip_list(edit_script_list **left)
{  
   edit_script_list  *ep, *ahead, *behind;
       
   ahead = *left;
   ep = NULL;
   while (ahead!=NULL) {
          behind = ep; ep = ahead;
          ahead = ahead->next_script;
          ep->next_script = behind;
  }  
  *left = ep;
}    


/* operates on a list sorted in increasing order of exon coordinates */
static void compact_list(Exon **Lblock, Exon **Rblock)
{
     Exon *tmp_block=*Lblock, *tmp_block1;
     int diff;

     while ((tmp_block!=NULL) && 
            ((tmp_block1=tmp_block->next_exon)!=NULL) &&
            tmp_block1->to1) {
        if ((abs((tmp_block1->from2-tmp_block1->from1) -
                (tmp_block->to2-tmp_block->to1))<=W) &&
                ((diff=tmp_block1->from2-tmp_block->to2-1)<=MAX_INTERNAL_GAP)) {
            /* merge blocks */
            tmp_block->to1 = tmp_block1->to1;
            tmp_block->to2 = tmp_block1->to2;
            tmp_block->length = tmp_block->to2-tmp_block->from2+1;
            tmp_block->edist += tmp_block1->edist;
            tmp_block->edist -= P*diff;
            tmp_block->next_exon = tmp_block1->next_exon;
    
            free(tmp_block1);
        } else
            tmp_block = tmp_block1;
     }
     /* reset right end of the list */
     *Rblock = tmp_block;
}

/* ------------------  memory management routines --------------- */

void link_to_data_list(Pointer data, ValNodePtr *head, ValNodePtr *prev)
{
     ValNodePtr curr;

     curr = (ValNodePtr)ckalloc(sizeof(struct ValNode));
     curr->data = data;
     curr->next = NULL;

     if(*prev == NULL)
        *head = curr;
     else
        (*prev)->next = curr;
     *prev = curr;
}

void   ValNodeFreeData(ValNodePtr data_list)
{
       ValNodePtr   tmp_node;

       while ((tmp_node=data_list)!=NULL) {
          free(tmp_node->data);
          data_list = data_list->next;
          free(tmp_node);
       }
}


int  good_ratio(int length)
{
     if (length<=W/2) return 2;
     else if (length<2*W) return rs.cutoff;
     else return (int)(.75*P*length+1);
}

static int extend_bw(uchar *s1, uchar *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)
{
  int     col,                    /* column number */
          row,                    /* row number */
          max_d,                  /* bound on the length of the edit script
 */
          d,                      /* current compressed distance */
          k,                      /* current diagonal */
          DELTA,                  /* n-m  */
          ORIGIN,
          lower,
          upper;
  int     *last_d, *temp_d;       /* column containing the last p */
  int     *min_row, *min_diag;    /* min (b)/ max (f) row (and diagonal) */
                                  /* reached for cost d=0, ... m.  */
  DELTA = n-m;
  max_d = m+1;

  ORIGIN = m;
  for (row=m, col=n; row>0 && col>0 && (s1[row-1]==s2[col-1]); row--,col--)
        /*LINTED empty loop body*/; 

  if ((row == 0) || (col == 0)) {
       *line1 = row+offset1;
       *line2 = col+offset2;

       return 0;
  }

  last_d = (int *)ckalloc((m+n+1)*sizeof(int));
  temp_d = (int *)ckalloc((m+n+1)*sizeof(int));

  for (k=0; k<=m+n; ++k) last_d[k]=m+1;
  last_d[ORIGIN+DELTA] = row;

  lower = ORIGIN + DELTA - 1;
  upper = ORIGIN + DELTA + 1;

  min_row = (int *)ckalloc((m+1)*sizeof(int));
  min_diag = (int *)ckalloc((m+1)*sizeof(int));

  for (d=1; d<=m; d++)
       min_row[d] = m+1;
 
  min_row[0] = last_d[ORIGIN+DELTA];
  min_diag[0] = ORIGIN + DELTA;

  d = 0;
  while ((++d<=max_d) && 
         ((d-1<=good_ratio(m-min_row[d-1])) ||
          ((d>=2) && (d-2<=good_ratio(m-min_row[d-2]))))) {

          /* for each relevant diagonal ... */
          for (k = lower; k <= upper; k++) {

               /* find a d on diagonal k */
               if (k==-d+DELTA+ORIGIN) {
                        /* move down from the last d-1 on diagonal k+1 */
                        row = last_d[k+1];
                        /* op = INSERT; */

               } else if (k==d+DELTA+ORIGIN) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1]-1;
                        /* op = DELETE; */

               } else if ((last_d[k]-1<=last_d[k+1]) &&
                          (last_d[k]-1<=last_d[k-1]-1)) {
                        /* substitution */
                        row = last_d[k]-1;
                        /* op = SUBSTITUTE; */

               } else if ((last_d[k-1]-1<=last_d[k+1]) &&
                          (last_d[k-1]-1<=last_d[k]-1)) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1]-1;
                        /* op = DELETE; */

               } else  {
                        /* move left from the last d-1 on diagonal k+1 */
                        row = last_d[k+1];
                        /* op = INSERT; */

               }

               /* code common to the three cases */
               /* slide down the diagonal */

               col = row+k-ORIGIN;

               while ((row > 0) && (col > 0) && (s1[row-1]==s2[col-1]))
                  { row--; col--; }

               temp_d[k] = row;

               if ((row == 0) && (col == 0)) {
                   /* hit southeast corner; have the answer */

                   free(last_d); free(temp_d);
                   free(min_row); free(min_diag);

                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }
               if (row == 0) {
                   /* hit first row; don't look further */

                   free(last_d); free(temp_d);
                   free(min_row); free(min_diag);

                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }

               if (col == 0) {
                   /* hit last column; don't look further */
                   free(last_d); free(temp_d);
                   free(min_row); free(min_diag);

                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }
          }

          min_row[d] = last_d[ORIGIN+DELTA];
          min_diag[d] = ORIGIN+DELTA;
          for (k=lower; k<=upper; ++k)
                 if (temp_d[k]<min_row[d]) {
                     min_row[d] = temp_d[k];
                     min_diag[d] = k;
                 }

          for (k=lower; k<=upper; k++) {
               last_d[k] = temp_d[k];
          }

          --lower;
          ++upper;
     }

    /* report here the previous maximal match, stored in min_diag and min_row */
     while ((d>0) && (min_row[d-1]-min_row[d]<3))
        d--;

     *line1 = min_row[d]+offset1;
     *line2 = min_row[d]+min_diag[d]-ORIGIN+offset2;

     free(min_row);
     free(min_diag);
     free(last_d);
     free(temp_d);

     return d;
}


static int extend_fw(uchar *s1, uchar *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)
{
  int     col,                    /* column number */
          row,                    /* row number */
          max_d,                  /* bound on the length of the edit script
 */
          d,                      /* current compressed distance */
          k,                      /* current diagonal */
          ORIGIN,
          lower,
          upper;
  int     *last_d, *temp_d;       /* column containing the last p */
  int     *max_row, *max_diag;    /* min (b)/ max (f) row (and diagonal) */
                                  /* reached for cost d=0, ... m.  */
  max_d = m+1;

  ORIGIN = m;
  for (row=0, col=0; col<n && row<m && (s1[row]==s2[col]); row++, col++)
        /*LINTED empty loop body*/; 

  if (row == m) {
       *line1 = row+offset1;
       *line2 = col+offset2;

       return 0;
  }
  if (col == n) {
       *line1 = row+offset1;
       *line2 = col+offset2;

       return 0;
  }

  last_d = (int *)ckalloc((m+n+1)*sizeof(int));
  temp_d = (int *)ckalloc((m+n+1)*sizeof(int));

  for (k=0; k<=m+n; ++k) last_d[k]=-1;
  last_d[ORIGIN] = row;

  lower = ORIGIN - 1;
  upper = ORIGIN + 1;

  max_row = (int *)ckalloc((m+1)*sizeof(int));
  max_diag = (int *)ckalloc((m+1)*sizeof(int));

  for (d=1; d<=m; d++)
       max_row[d] = -1;
 
  max_row[0] = last_d[ORIGIN];
  max_diag[0] = ORIGIN;


  d = 0;
  while ((++d<=max_d) && 
         ((d-1<=good_ratio(max_row[d-1])) || 
          ((d>=2) && (d-2<=good_ratio(max_row[d-2]))))) {

          /* for each relevant diagonal ... */
          for (k = lower; k <= upper; k++) {

               /* find a d on diagonal k */
               if (k==-d+ORIGIN) {

                        /* move down from the last d-1 on diagonal k+1 */
                        row = last_d[k+1]+1;
                        /* op = DELETE; */
               } else if (k==d+ORIGIN) {

                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1];
                        /* op = INSERT; */
               } else if ((last_d[k]>=last_d[k+1]) &&
                          (last_d[k]+1>=last_d[k-1])) {

                        /* substitution */
                        row = last_d[k]+1;
                        /* op = SUBSTITUTE; */
               } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                          (last_d[k+1]>=last_d[k])) {

                        /* move down from the last d-1 on diagonal k+1 */
                        row = last_d[k+1]+1;
                        /* op = DELETE; */
               } else {

                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1];
                        /* op = INSERT; */
               }

               /* code common to the three cases */
               /* slide down the diagonal */

               col = row+k-ORIGIN;

               if (row>=0)
               while ((row < m) && (col < n) && (s1[row]==s2[col]))
                         { row++; col++; }

               temp_d[k] = row;

               if ((row == m) && (col == n)) {
                   /* hit southeast corner; have the answer */
                   free(last_d); free(temp_d);
                   free(max_row); free(max_diag);
                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }
               if (row == m) {
                   /* hit last row; don't look further */
                   free(temp_d); free(last_d);
                   free(max_row); free(max_diag);

                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }

               if (col == n) {
                   /* hit last column; don't look further */
                   free(temp_d); free(last_d);
                   free(max_row); free(max_diag);

                   *line1 = row+offset1;
                   *line2 = col+offset2;

                   return d;
               }
          }
          max_row[d] = last_d[ORIGIN];
          max_diag[d] = ORIGIN;
          for (k=lower; k<=upper; ++k)
                 if (temp_d[k]>max_row[d]) {
                     max_row[d] = temp_d[k];
                     max_diag[d] = k;
                 }

          for (k=lower; k<=upper; k++) {
               last_d[k] = temp_d[k];
          }

          --lower;
          ++upper;
     }

    /* report here the previous maximal match, stored in max_diag and max_row */

     while ((d>0) && (max_row[d]-max_row[d-1]<3))
        d--;

     *line1 = max_row[d]+offset1;
     *line2 = max_row[d]+max_diag[d]-ORIGIN+offset2;

     free(max_row);
     free(max_diag);
     free(last_d);
     free(temp_d);

     return d;

/*
     if ((d>2) && (max_row[d-1]-max_row[d-2]<3)) {
          *line1 = max_row[d-2]+offset1;
          *line2 = max_row[d-2]+max_diag[d-2]-ORIGIN+offset2;
 
          free(max_row); free(max_diag);
          free(last_d); free(temp_d);
 
          return d-2;
     }
 
     *line1 = max_row[d-1]+offset1;
     *line2 = max_row[d-1]+max_diag[d-1]-ORIGIN+offset2;

     free(max_row);
     free(max_diag);
     free(last_d);
     free(temp_d);

     return d-1;
*/
}

static  void merge(Exon **t0, Exon **t1) 
{   
    Exon  *tmp0, *tmp1;
    int    diff;

    if ((*t0) && !(*t0)->to1)
         tmp0 = (*t0)->next_exon;
    else    
         tmp0 = *t0;

    while (tmp0 && (tmp0!=*t1)) {
       tmp1 = tmp0->next_exon;
       assert(tmp1!=NULL);
       if (tmp1 && tmp1->to1 && tmp0->to1 && 
           (abs((tmp1->from2-tmp1->from1)-(tmp0->to2-tmp0->to1))<=W) &&
           ((diff=tmp1->from2-tmp0->to2-1<=W))) {
         
          /* merge blocks tmp0 and tmp1 */
          tmp0->from1 = min(tmp0->from1, tmp1->from1);
          tmp0->from2 = min(tmp0->from2, tmp1->from2);
          tmp0->to1 = max(tmp1->to1, tmp0->to1);
          tmp0->to2 = max(tmp1->to2, tmp0->to2);
          tmp0->length = tmp0->to2-tmp0->from2+1;
          tmp0->flag = tmp1->flag; 
          tmp0->edist += tmp1->edist; 
          tmp0->edist -= P*diff;
          if (tmp1==*t1) {
          /*  tmp0->flag = (*t1)->flag; */
              *t1 = tmp0;
          }
          tmp0->next_exon = tmp1->next_exon;  

          free(tmp1);
       } else 
          tmp0 = tmp0->next_exon;
    }
}

void free_align(edit_script_list *aligns)
{
     edit_script_list *head;

     head = aligns;

     while ((head=aligns)!=NULL) {
        aligns = aligns->next_script;
        Free_script(head->script);
        free(head); 
     }
}

     
void free_list(Exon *left)
{
     Exon *tmp_block;
      
     while ((tmp_block=left)!=NULL) {
          left = left->next_exon;
          free(tmp_block);
     }
}    

void free_table(void)
{    
     register struct hash_node *hptr, *tptr;
     register int    hval;

     free(pnext_pos);
     for (hval=0; hval<HASH_SIZE+1; hval++) {
          hptr = phashtab[hval];
          while (hptr) {
             tptr = hptr;
             hptr = hptr->link;
             free(tptr);
          } 
     }     
}  

static Exon *bmatch (uchar *s1, uchar *s2, int len1, int len2, int offset1, int offset2)
{
       int  i, j, i1, score;
       Exon *new=NULL;
 
       for (i1=i=len1-3; i>=len2-3; i--, i1=i) {
         for (j=len2-3; j>=2; j--, i1--)
           if (*(s1+i1)!=*(s2+j))
               break;
          
         if (j<2) {
            /* exact match for CDS found; check signals */
            score = 0;
            if (*(s1+(i1--))==*(s2+(j--))) score++;
            if (*(s1+(i1--))==*(s2+(j--))) score++;
            if (*(s1+i1+len2-1)==*(s2+j+len2-1)) score++;
            if (*(s1+i1+len2)==*(s2+j+len2)) score++;
            if (score>=3) {
                new = new_exon(i1+3+offset1, offset2, i1+3+offset1+len2-5,
                               offset2+len2-5, len2-4, 0, 0, NULL);
                new->ori = (G_score >= abs(C_score)) ? 'G' : 'C';
    
                return new;
            }
         }
       }
       return NULL;
}

static Exon *fmatch (uchar *s1, uchar *s2, int len1, int len2, int offset1, int offset2)
{
       int  i, j, i1, score;
       Exon *new=NULL;
 
       for (i1=i=2; i<len1-len2+3; i++, i1=i) {
         for (j=2; j<len2-2; j++, i1++)
           if (*(s1+i1)!=*(s2+j))
               break;

         if (j>=len2-2) {
            /* exact match found for internal part, look for signals */
            score = 0; 
            if (*(s1+(i1++))==*(s2+(j++))) score++;
            if (*(s1+(i1++))==*(s2+(j++))) score++;
            if (*(s1+i1-len2)==*s2) score++;
            if (*(s1+i1-len2+1)==*(s2+1)) score++;
            if (score>=3) {
                new = new_exon(i+offset1,offset2,i1+offset1-2,offset2+len2-5,
                               len2-4,0,0,NULL);
                new->ori = (G_score >= abs(C_score)) ? 'G' : 'C';

                return new;
            }
         }  
       }    
       return NULL;
}

#ifdef DEBUG
static void debug_print_exons (Exon *lblock, const char *label)
{
   Exon *tmp_block = lblock;

   (void)fprintf(stderr,"\n====================%s:\n\n", label);
   while (tmp_block) {
      (void)fprintf(stderr," [ %d, %d, %d, %d, l: %d ]\n ", tmp_block->from1, 
                    tmp_block->from2, tmp_block->to1,
                    tmp_block->to2, tmp_block->length);
      tmp_block = tmp_block->next_exon;
  }  
}
#endif

/* -------------------- to be added to psublast ---------------------- */

void seq_toupper(uchar *seq, int len, char *filename)
{
     int i=0, flag = 0;
     uchar *s=seq;

     for (; *s && (i<min(100,len)); i++, s++) 
       if (islower(*s)) { flag = 1; break; }

     if (flag) { 
       fprintf(stderr, "Warning: lowercase letter in %s data. Convert sequence to uppercase.\n", filename);
       for (s=seq, i=0; *s && (i<len); s++, i++)
            *s = toupper(*s);
     }
}

static bool get_sync_flag(Exon *lblock, Exon *rblock, int w)
{
     int numx=0, e2;
     Exon *t;

     if (((t=lblock->next_exon)==NULL) || !t->to1) return FALSE; 
     numx++; e2 = t->to2; 

     while (((t=t->next_exon)!=NULL) && t->to1) {
        ++numx;
        if ((t->from2-e2>1) || 
            (t!=rblock && ((t->to2-t->from2+1<2*w+2) ||
                           (t->to1-t->from1+1<2*w+2)))) return FALSE;
        e2 = t->to2;
     }
     return ((numx<3) ? FALSE:TRUE);
}
 

static void sync_slide_intron(int in_w, Exon **lblock, uchar *seq1, uchar *seq2)
{
    Exon *t0=NULL, *t1=NULL, *head = *lblock;
    splice_t *g=NULL, *c=NULL, *cell=NULL;
    splice_t *Glist[500], *Clist[500];
    int Gscore=0, Cscore=0;
    char oris[500];
    int w1, w2, ni, i, numC, numG;

    memset(Glist, 0, 200*sizeof(splice_t *));
    memset(Clist, 0, 200*sizeof(splice_t *));

    ni = 0; numG = numC = 0;

    /* assume forward orientation */
    t0 = head->next_exon;
    while (t0 && (t1=t0->next_exon) && t1->to1) {
       g = c = NULL;
       if (t1->from2-t0->to2-1==0) {
         if (!strncmp((char *)(seq1+t0->to1),"GT",2) &&
             !strncmp((char *)(seq1+t1->from1-3),"AG",2)) {
              g = new_splice('G',t0->to1,t1->from1,t0->to2,t1->from2,-1,NULL);
              t0->ori = 'G';
              oris[ni] = 'G'; 
              numG++;
         } else if (!strncmp((char *)(seq1+t0->to1),"CT",2) &&
                  !strncmp((char *)(seq1+t1->from1-3),"AC",2)) {
              c = new_splice('C',t0->to1,t1->from1,t0->to2,t1->from2,-1,NULL);
              t0->ori = 'C';  
              oris[ni] = 'C';
              numC++;
         } else {
              w1 = min(in_w, min(t0->length-1, t0->to1-t0->from1));
              w2 = min(in_w, min(t1->length-1, t1->to1-t1->from1));
              splice(seq1, t0->to1-w1, t0->to1+w1, t1->from1-w2, t1->from1+w2,
                     seq2, t0->to2-w1, t1->from2+w2, &g, &c, BOTH);

              Gscore += g->score; Cscore += c->score;
              cell = NULL; oris[ni] = '*';
              if (g->score>c->score) {
                  numG++; cell = g; oris[ni] = 'G';
              } else if (c->score>g->score) {
                  numC++; cell = c; oris[ni] = 'C';
              } else if (c->score==g->score) {
                  numG++; numC++; cell = g;  oris[ni] = 'G';
              }
              t0->ori = oris[ni];
              t0->to1 = cell->xs; t0->to2 = cell->ys;
              t1->from1 = cell->xe; t1->from2 = cell->ye;
              t0->length = t0->to2-t0->from2+1;
              t1->length = t1->to2-t1->from2+1;
         }
         Clist[ni] = c; Glist[ni] = g;
       } else {
         t0->ori = 'E'; oris[ni] = 'E';
       }
       ni++;
       t0 = t1;
    }         
              
/*
    if (numG==ni) {
        for (i=0; i<ni; i++) {
             if (Glist[i]) free(Glist[i]);
             if (Clist[i]) free(Clist[i]);
        }     
        return;
    }         
*/
    if ((numG==1) && (numC==1) && 
        (!Glist[0] || !Clist[0] || !Glist[1] || !Clist[1])) goto free_all;
              
    if (numG>=numC) {
        /* revisit all previous assignments that are inconsistent */
        for (i=0, t0=head->next_exon; i<ni; i++, t0=t1) {
          t1 = t0->next_exon;
          switch (oris[i]) {
             case 'G': break;
             case 'C': if (Glist[i]==NULL) { 
                          /* compute the values for C */ 
                          w1 = min(in_w, min(t0->length-1, t0->to1-t0->from1));
                          w2 = min(in_w, min(t1->length-1, t1->to1-t1->from1));
                          splice(seq1, t0->to1-w1, t0->to1+w1,
                                 t1->from1-w2, t1->from1+w2, seq2,
                                 t0->to2-w1, t1->from2+w2, &g, &c, FWD);
                       } else g = Glist[i];
              
                       t0->ori = 'G';
                       t0->to1 = g->xs; t0->to2 = g->ys;
                       t1->from1 = g->xe; t1->from2 = g->ye;
                       t0->length = t0->to2-t0->from2+1;
                       t1->length = t1->to2-t1->from2+1;

                       break;
             case 'E': break;
             default : fatal("sim4b1.c: intron orientation not initialized.");
          }
          if (oris[i]!='E') wobble(&t0,&t1,"GT","AG",seq1);
        }
    } else {
        /* analyze all assignments for consistency */
        for (i=0, t0=head->next_exon; i<ni; i++, t0=t1) {
          t1 = t0->next_exon;
          switch (oris[i]) {
             case 'C': break;
             case 'G': if (Clist[i]==NULL) {
                          /* compute the values for C */
                          w1 = min(in_w, min(t0->length-1, t0->to1-t0->from1));
                          w2 = min(in_w, min(t1->length-1, t1->to1-t1->from1));
                          splice(seq1, t0->to1-w1, t0->to1+w1,
                                 t1->from1-w2, t1->from1+w2,
                                 seq2, t0->to2-w1, t1->from2+w2, &g, &c, BWD);
                       } else c = Clist[i];
          
                       t0->ori = 'C';
                       t0->to1 = c->xs; t0->to2 = c->ys;
                       t1->from1 = c->xe; t1->from2 = c->ye;
                       t0->length = t0->to2-t0->from2+1; 
                       t1->length = t1->to2-t1->from2+1;
                       break;
             case 'E': break;
             default : fatal("sim4b1.c: intron orientation not initialized.");
          }
          if (oris[i]!='E') wobble(&t0,&t1,"CT","AC",seq1);
        }
    }

    /* now free all memory allocated */
    free_all:
    for (i=0; i<ni; i++) {
        if (Glist[i]) free(Glist[i]);
        if (Clist[i]) free(Clist[i]);
    }
    return;
}

static void wobble(Exon **t0, Exon **t1, const char *donor, const char *acceptor, uchar *seq1)
{
    uchar *s = seq1+(*t0)->to1;  /* first nt of donor */
    uchar *q = seq1+(*t1)->from1-3;  /* first nt of acceptor */

    if (!strncmp((char *)(s), donor, 2)) {
       /* match in place */
       if (!strncmp((char *)(q), acceptor, 2)) {
          return; 
       } else if (!strncmp((char *)(q-1), acceptor, 2)) {
          (*t1)->from1--; return;
       } else if (!strncmp((char *)(q+1), acceptor, 2)) {
          (*t1)->from1++; return;
       }

    } else if (!strncmp((char *)(s-1), donor, 2)) {
       /* match is 1 off to the left */
       if (!strncmp((char *)(q), acceptor, 2)) {   
          (*t0)->to1--; return;
       } else if (!strncmp((char *)(q-1), acceptor, 2)) {
          (*t0)->to1--; (*t1)->from1--; 
          (*t0)->to2--; (*t1)->from2--;
          (*t0)->length++; (*t1)->length--;
          return;
       } else if (!strncmp((char *)(q+1), acceptor, 2)) {
          (*t0)->to1--; (*t1)->from1++; return;
       }
       
    } else if (!strncmp((char *)(s+1), donor, 2)) {
       /* match is 1 off to the right */
       if (!strncmp((char *)(q), acceptor, 2)) {
          (*t0)->to1++; return;
       } else if (!strncmp((char *)(q-1), acceptor, 2)) {
          (*t0)->to1++; (*t1)->from1--; return;
       } else if (!strncmp((char *)(q+1), acceptor, 2)) {
          (*t0)->to1++; (*t1)->from1++; 
          (*t0)->to2++; (*t1)->from2++;
          (*t0)->length--; (*t1)->length++;
          return;
       }
    } else if (!strncmp((char *)(q-1), acceptor, 2)) {
       /* match is 1 off to the left */
       (*t1)->from1--; return;

    } else if (!strncmp((char *)(q+1), acceptor, 2)) {
       /* match is 1 off to the right */
          (*t1)->from1++; return;
    }

    return;
}

static void slide_intron(int in_w, Exon **lblock, uchar *seq1, uchar *seq2)
{ 
    Exon *t0, *t1, *head = *lblock;       
    splice_t *g, *c, *cell;
    char type;
    int w1, w2;

    t0 = head->next_exon;
    while (t0 && (t1=t0->next_exon) && t1->to1) {
       g = c = NULL;
       if (t1->from2-t0->to2-1==0) {
         if (!strncmp((char *)(seq1+t0->to1),"GT",2) && 
             !strncmp((char *)(seq1+t1->from1-3),"AG",2)) 
                   t0->ori = 'G'; 
         else if (!strncmp((char *)(seq1+t0->to1),"CT",2) && 
                  !strncmp((char *)(seq1+t1->from1-3),"AC",2))
                   t0->ori = 'C';
         else {      
           int gtag=0, ctac=0;
           uchar *s;
   
           w1 = min(in_w, min(t0->length-1, t0->to1-t0->from1));
           w2 = min(in_w, min(t1->length-1, t1->to1-t1->from1));
           splice(seq1, t0->to1-w1, t0->to1+w1, t1->from1-w2, t1->from1+w2,
                  seq2, t0->to2-w1, t1->from2+w2, &g, &c, BOTH);
           if (g->score>c->score) { cell = g; type = 'G'; }
           else if (c->score>g->score) { cell = c; type = 'C'; }
           else { cell = g; type = 'G'; }
 
           t0->to1 = cell->xs; t0->to2 = cell->ys;  
           t1->from1 = cell->xe; t1->from2 = cell->ye;   
           t0->length = t0->to2-t0->from2+1;   
           t1->length = t1->to2-t1->from2+1;   
           
           wobble(&t0,&t1,(type=='G')? "GT":"CT",(type=='G')? "AG":"AC",seq1);

           free(g); free(c);

           /* determine the type, based on the # matches w/ GT-AG (CT-AC) */
           s = seq1+t0->to1;
           if (*s=='G') gtag++; else if (*s=='C') ctac++; 
           ++s;
           if (*s=='T') { gtag++; ctac++;} 
           s = seq1+t1->from1-3; 
           if (*s=='A') { gtag++; ctac++; }
           ++s;
           if (*s=='G') gtag++; else if (*s=='C') ctac++;
           if (gtag>ctac) type = 'G';
           else if (ctac>gtag) type = 'C';
           else type = 'N';

           t0->ori = type;         
         }
       } else
           t0->ori = 'E';
       t0 = t1;
    }      
    
}

#ifdef AUXUTILS
static void remove_polyA_tails(Exon *lblock, uchar *seq1, uchar *seq2, int len2)
{
    Exon *t, *prev;
    uchar *s, *q;
    int xcut, diff, u, tmp, I, J, first = 1;

    t = lblock->next_exon; prev = lblock;

    s = seq2; q = seq2+len2-1; while (s<=q && *s=='T') s++;
    if (s>=seq2+t->from2-1) {
        while (t && t->to1) {
           s = seq2+t->from2-1; q = seq2+t->to2-1;
           if (first && strncmp((char *)s,"TTTTT",5)) break;
           first = 0;
           while ((s<=q) && (*s=='T')) s++;
           diff = t->to2-(s-seq2); u = min(diff,12); 
           xcut = t->to2-t->from2+1-diff; 
           if (diff>6) {
               tmp = (int)((1+P)*u); /* was diff */
/*
               EXTEND_BW(s,seq1+t->to1-(diff-u)-tmp,u,tmp,
                         s-seq2,t->to1-(diff-u)-tmp,&I,&J);
*/
               EXTEND_BW(seq2+t->from2+xcut-1,
                         seq1+t->from1+(int)((1-P)*xcut)-1,
                         u,(int)(P*xcut+(1+P)*u)+1, /* 1 is for looser margin */
                         t->from2+xcut-1,t->from1+(int)((1-P)*xcut)-1,
                         &I,&J);
               t->from2 = I+1; t->from1 = J+1;
               t->length = t->to2-t->from2+1;
               break;
           } else if (diff>0) {
               prev->next_exon = t->next_exon; free(t); t = prev->next_exon;
    
               tmp = (int)((1+P)*diff);
               EXTEND_BW(s,seq1+t->from1-tmp,diff,tmp,
                         s-seq2,t->from1-tmp,&I,&J);
               t->from2 = I+1; t->from1 = J+1; t->length = t->to2-t->from2+1;
               break;
           } else {
               /* remove entire exon and repeat the process */
               prev->next_exon = t->next_exon; free(t); t = prev->next_exon;
               continue;
           }
        }
    }

    while (t && t->next_exon && t->next_exon->to1) { 
       prev = t; t = t->next_exon; 
    }

    first = 1;
    s = seq2; q = seq2+len2-1; while (q>=s && *q=='A') q--;
    if (t && t->to1 && q<seq2+t->to2-1) {
        while (t && t->to1) {
           s = seq2+t->to2-1; q = seq2+t->from2-1;
           if (first && strncmp((char *)(s-4), "AAAAA", 5)) break;
           first = 0;
           while ((s>=q) && (*s=='A')) s--;
           diff = (int)(s-seq2+1)-t->from2+1; u = min(diff, 12);
           xcut = t->to2-t->from2+1-diff;
           if (diff>6) {
/*
               EXTEND_FW(seq2+t->from2+(diff-u)-1, seq1+t->from1+(diff-u)-1, 
                         u, (int)((1+P)*u), t->from2+(diff-u)-1, 
                         t->from1+(diff-u)-1, &I, &J);

*/
               EXTEND_FW(s-u+1, seq1+t->to1-(int)((1+P)*xcut)-(int)((1+P)*u), 
                         u,(int)((1+P)*u+P*xcut)+1, /*+1 is for looser margin */
                         t->from2+(diff-u)-1,
                         t->to1-(int)((1+P)*xcut)-(int)((1+P)*u), &I, &J);

               t->to1 = J; t->to2 = I; t->length = t->to2-t->from2+1;
               break;
           } else if (diff>0) {
               if (prev==NULL) prev = find_previous(lblock, t);
               assert(prev!=NULL);
               prev->next_exon = t->next_exon; free(t); 
               t = prev->next_exon; prev = NULL; 
               EXTEND_FW(seq2+t->to2, seq1+t->to1, diff, (int)((1+P)*diff),
                         t->to2, t->to1, &I, &J);
               t->to1 = J; t->to2 = I; t->length = t->to2-t->from2+1;
               break;
           } else {
               if (prev==NULL) prev = find_previous(lblock, t);
               assert(prev!=NULL);
               prev->next_exon = t->next_exon; free(t); 
               t = prev->next_exon; prev = NULL;
               continue;
           }
       }
    }

    return;
}
#endif

Exon *find_previous(Exon *lblock, Exon *tgt)
{
    Exon *t=lblock;

    while (t && (t->next_exon!=tgt)) t = t->next_exon; 
    if (t==NULL) 
        fatal("sim4b1.c: Corrupted exon list: could not find previous.");

    return t;
}    

static bool get_match_quality(Exon *lblock, Exon *rblock, sim4_stats_t *st, int N)
{
  int  tcov;
  bool good_match;
  Exon *t;
  
  good_match = TRUE; st->icoverage = 0;
  t = lblock->next_exon;
  while (t->to1) {
     st->icoverage += t->to2-t->from2+1;
     if (100*t->edist>=5*(t->to2-t->from2+1)) {
         good_match = FALSE; break;
     }
     t = t->next_exon; 
  }
  tcov = rblock->to2-lblock->next_exon->from2+1;
  if (lblock->next_exon->from2>=.5*N &&
      tcov>=.8*(N-lblock->next_exon->from2) &&
      st->icoverage>=max(.95*tcov,100))
         ;
  else if (rblock->to2<=.5*N && tcov>=.8*rblock->to2 &&
           st->icoverage>=max(.95*tcov,100))
         ;
  else if ((tcov<.8*N) || (st->icoverage<.9*tcov))
                 good_match = FALSE;
  
  return good_match; 
}
