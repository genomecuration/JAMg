#include "libc.h"
#include "sim4.h"
#include "psublast.h"
#include "splice.h"

#ifndef __lint
/*@unused@*/
static const char rcsid[] =
"$Id: splice.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif


static int spl_encode[NACHARS];
static int encodeInit;

signal_t gt = {{0, 0, 0, 2, 0},{0, 0, 0, 2, 0},{2, 3, 2, 5, 2},{0, 0, 0, 2, 0},{0, 0, 0, 2, 0}};
signal_t ct = {{0, 0, 0, 2, 0},{2, 2, 2, 5, 2},{0, 0, 0, 2, 0},{0, 0, 0, 2, 0},{0, 0, 0, 2, 0}};
signal_t ag = {{2, 2, 5, 2, 2},{0, 0, 2, 0, 0},{0, 0, 2, 0, 0},{0, 0, 2, 0, 0},{0, 0, 2, 0, 0}};
signal_t ac = {{2, 5, 2, 2, 2},{0, 2, 0, 0, 0},{0, 3, 0, 0, 0},{0, 2, 0, 0, 0},{0, 2, 0, 0, 0}};


#ifdef DEBUG
static void print_splice(splice_t *);
#endif

static void splice_donor(uchar *xseq, uchar *yseq, int M, int N, int *gt_score,
                         int *ct_score, int **max_Gf, int **max_Cf,
                         int **start_Gi, int **start_Ci);
static void splice_donor_uni(uchar *xseq, uchar *yseq, int M, int N,
                         int *It_score, int **max_IF, int **end_Ii);
static void splice_acceptor(uchar *xseq, uchar *yseq, int M, int N,
                         int *ag_score, int *ac_score, int **max_Gb,  
                         int **max_Cb, int **end_Gi, int **end_Ci);
static void splice_acceptor_uni(uchar *xseq, uchar *yseq, int M, int N,
                         int *aI_score, int **max_Ib, int **end_Ii);
static int  stepct(int n);


void splice(uchar *in_seqx, int ls, int us, int le, int ue,
            uchar *in_seqy, int ys, int ye, 
            splice_t **gcell, splice_t **ccell, int ori)
{
   int p, q, *gtscore=NULL, *ctscore=NULL, *agscore=NULL, *acscore=NULL;
   int i, tmp;
   int maxCscore, maxGscore, Gxs, Gxe, Gy, Cxs, Cxe, Cy, keep_Ci, keep_Gi;
   int *max_Cf=NULL, *max_Gf=NULL, *max_Cb=NULL, *max_Gb=NULL;
   int *start_Gi=NULL, *start_Ci=NULL, *end_Gi=NULL, *end_Ci=NULL;
   uchar *s;
   
   if (!encodeInit) {
       for (i=0; i<NACHARS; spl_encode[i++]=4);
     
       spl_encode['A'] = 0; spl_encode['C'] = 1;
       spl_encode['G'] = 2; spl_encode['T'] = 3;

       encodeInit = 1;
   }
       
   if (ori==FWD || ori==BOTH) {
       gtscore = (int *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(int));
       agscore = gtscore+(us-ls+2);
       for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++) 
            gtscore[p] = gt[spl_encode[*s]][spl_encode[*(s+1)]];
       for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
            agscore[q] = ag[spl_encode[*(s-1)]][spl_encode[*s]]; 
   } 
   if (ori==BWD || ori==BOTH) {
       ctscore = (int *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(int));
       acscore = ctscore+(us-ls+2);
       for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++)  
            ctscore[p] = ct[spl_encode[*s]][spl_encode[*(s+1)]];  
       for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
            acscore[q] = ac[spl_encode[*(s-1)]][spl_encode[*s]];
   } 

   if (ori==FWD) {
      splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                       gtscore, &max_Gf, &start_Gi);
      splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                       agscore, &max_Gb, &end_Gi);
      free(gtscore); /* free(agscore); */

   } else if (ori==BWD) {
      splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                       ctscore, &max_Cf, &start_Ci);              
      splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                       acscore, &max_Cb, &end_Ci);
      free(ctscore); /* free(acscore); */

   } else {
      splice_donor(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1, 
                   gtscore, ctscore, &max_Gf, &max_Cf, &start_Gi, &start_Ci);

      splice_acceptor(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1, 
                      agscore, acscore, &max_Gb, &max_Cb, &end_Gi, &end_Ci);

      free(gtscore); /* free(agscore);     */
      free(ctscore); /* free(acscore); */
   }

   maxCscore = -999999; maxGscore = -999999;
   Gxs = Gxe = Gy = Cxs = Cxe = Cy = -1;
   if (ori==FWD || ori==BOTH) {
      for (i=0; i<=ye-ys+1; i++) {
        if ((tmp=max_Gf[i]+max_Gb[i])>maxGscore) {
            maxGscore = tmp;
            /* save (i, start_Gi[i], end_Gi[i]); */
            Gxs = ls+start_Gi[i]-1; Gxe = le+end_Gi[i]-1; Gy = ys+i-1;
            keep_Gi = i;
        }  
      }
      free(max_Gf); free(max_Gb); /* free(start_Gi); free(end_Gi); */
   }
   if (ori==BWD || ori==BOTH) {
      for (i=0; i<=ye-ys+1; i++) {
        if ((tmp=max_Cf[i]+max_Cb[i])>maxCscore) {
            maxCscore = tmp;
            /* save (i, start_Ci[i], end_Ci[i]); */
            Cxs = ls+start_Ci[i]-1; Cxe = le+end_Ci[i]-1; Cy = ys+i-1;
            keep_Ci = i;
        }
      }
      free(max_Cf); free(max_Cb); /* free(start_Ci); free(end_Ci); */
   }

   *gcell = new_splice('G', Gxs, Gxe, Gy, Gy+1, maxGscore, NULL);
   *ccell = new_splice('C', Cxs, Cxe, Cy, Cy+1, maxCscore, NULL);

#ifdef DEBUG
   print_splice(*gcell); print_splice(*ccell);
#endif

   return;
}

splice_t *new_splice(char c, int xs, int xe, int ys, int ye, int score, splice_t *next)
{
   splice_t *sp = ckalloc(sizeof(splice_t));
   
   sp->type = c; sp->xs = xs; sp->xe = xe; 
   sp->ys = ys; sp->ye = ye; sp->score = score;
   sp->next = next;

   return sp; 
}

#ifdef DEBUG
static void print_splice(splice_t *g)
{
    printf("Type: %c  sx: %d  se: %d  ys: %d  score: %d\n",
            g->type, g->xs, g->xe, g->ys, g->score);

}
#endif

static void splice_donor(uchar *xseq, uchar *yseq, int M, int N, int *gt_score,
                         int *ct_score, int **max_Gf, int **max_Cf, 
                         int **start_Gi, int **start_Ci)
{
   int *CCf, *mG, *mC, *sC, *sG, *X;
   int i, j, tmp, ss, ssx, cx, c;
   uchar *s, *t;

   CCf = (int *)ckalloc((M+1)*sizeof(int));
   X = (int *)ckalloc((M+1)*sizeof(int));
   mG = *max_Gf = (int *)ckalloc((2*N+2)*sizeof(int));
   sG = *start_Gi = mG+(N+1);
   mC = *max_Cf = (int *)ckalloc((2*N+2)*sizeof(int));
   sC = *start_Ci = mC+(N+1);

   t = yseq; X[0] = CCf[0] = 0;
   for (j=1; j<=M; j++) { CCf[j] = j; X[j] = 0; } 
   
   mG[0] = mC[0] = -999999;
   for (j=0; j<=M; j++) {
      if ((100*gt_score[j])>mG[0]) { mG[0] = 100*gt_score[j]; sG[0] = j; }
      if ((100*ct_score[j])>mC[0]) { mC[0] = 100*ct_score[j]; sC[0] = j; } 
   }
   
   for (i=1; i<=N; i++, t++) {
     s = xseq;
     ss = CCf[0]; ssx = X[0];
     c = ++CCf[0]; cx = X[0];
     for (j=1; j<=M; j++, s++) {
          tmp=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
          if (tmp==c+1);
          else if (tmp==CCf[j]+1) cx = X[j];
          else cx = ssx + (*t==*s);
          c = tmp; ss = CCf[j]; CCf[j] = c; ssx = X[j]; X[j] = cx;
     }
   
     /* compute max_Gf and max_Cf */
     mG[i] = mC[i] = -999999; 
     for (j=0; j<=M; j++) {
        assert(X[j]+CCf[j]!=0);
        tmp = (int)(stepct(j)*X[j]/(double)(X[j]+CCf[j])*100); 
        if ((tmp+100*gt_score[j])>mG[i]) {
            mG[i] = tmp+100*gt_score[j]; sG[i] = j;
        }
        if ((tmp+100*ct_score[j])>mC[i]) {
            mC[i] = tmp+100*ct_score[j]; sC[i] = j;
        }
     } 
   } 
   free(CCf); free(X); 
}  

static void splice_donor_uni(uchar *xseq, uchar *yseq, int M, int N,
                             int *It_score, int **max_If, int **start_Ii)
{
   int *CCf, *mI, *sI, *X;
   int i, j, tmp, ss, ssx, cx, c;
   uchar *s, *t;

   CCf = (int *)ckalloc((M+1)*sizeof(int));
   X = (int *)ckalloc((M+1)*sizeof(int));
   mI = *max_If = (int *)ckalloc((2*N+2)*sizeof(int));
   sI = *start_Ii = mI+(N+1);
  
   t = yseq; X[0] = CCf[0] = 0;
   for (j=1; j<=M; j++) { CCf[j] = j; X[j] = 0; }
 
   mI[0] = -999999;        
   for (j=0; j<=M; j++)                  
      if ((100*It_score[j])>mI[0]) { mI[0] = 100*It_score[j]; sI[0] = j; }      
 
   for (i=1; i<=N; i++, t++) {
     s = xseq;
     ss = CCf[0]; ssx = X[0];
     c = ++CCf[0]; cx = X[0];
     for (j=1; j<=M; j++, s++) {
          tmp=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
          if (tmp==c+1);
          else if (tmp==CCf[j]+1) cx = X[j];
          else cx = ssx + (*t==*s);
          c = tmp; ss = CCf[j]; CCf[j] = c; ssx = X[j]; X[j] = cx;
     }
  
     /* compute max_If */           
     mI[i] = -999999;        
     for (j=0; j<=M; j++) {
        assert(X[j]+CCf[j]!=0);
        tmp = (int)(stepct(j)*X[j]/(double)(X[j]+CCf[j])*100)+100*It_score[j];
        if (tmp>mI[i]) {
            mI[i] = tmp; sI[i] = j;
        }
     }
   }
   free(CCf); free(X);
}

     
static void splice_acceptor(uchar *xseq, uchar *yseq, int M, int N, 
                            int *ag_score, int *ac_score, int **max_Gb, 
                            int **max_Cb, int **end_Gi, int **end_Ci)
{
   int *CCb, *X, *mC, *mG, *eC, *eG;
   int tmp, i, j, ss, ssx, cx, c;
   uchar *t, *s;

   CCb = (int *)ckalloc((M+1)*sizeof(int));
   X = (int *)ckalloc((M+1)*sizeof(int));
   mG = *max_Gb = (int *)ckalloc((2*N+2)*sizeof(int));
   eG = *end_Gi = mG+(N+1);
   mC = *max_Cb = (int *)ckalloc((2*N+2)*sizeof(int));
   eC = *end_Ci = mC+(N+1);
   
   t = yseq+N-1; CCb[M] = X[M] = 0;
   for (j=M-1; j>=0; j--) { CCb[j] = M-j; X[j] = 0; }
   
   mG[N] = mC[N] = -999999;
   for (j=M; j>=0; j--) {
      if ((100*ag_score[j])>mG[N]) { mG[N] = 100*ag_score[j]; eG[N] = j+1; }
      if ((100*ac_score[j])>mC[N]) { mC[N] = 100*ac_score[j]; eC[N] = j+1; }
   }
   
   for (i=N-1; i>=0; i--, t--) {
     s = xseq+M-1; 
     ss = CCb[M]; ssx = X[M];
     c = ++CCb[M]; cx = X[M];
     for (j=M-1; j>=0; j--, s--) {
          tmp=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
          if (tmp==c+1) ; 
          else if (tmp==CCb[j]+1) cx = X[j];
          else cx = ssx + (*t==*s);
          c = tmp; ss = CCb[j]; CCb[j] = c; ssx = X[j]; X[j] = cx;
     }
     
     /* compute max_Gb and max_Cb */
     mG[i] = -999999; mC[i] = -999999;
     for (j=M; j>=0; j--) {
        assert(CCb[j]+X[j]!=0);
        tmp = (int)(stepct(M-j)*X[j]/(double)(CCb[j]+X[j])*100);
        if ((tmp+100*ag_score[j])>mG[i]) {
            mG[i] = tmp+100*ag_score[j]; eG[i] = j+1;
        }
        if ((tmp+100*ac_score[j])>mC[i]) {
            mC[i] = tmp+100*ac_score[j]; eC[i] = j+1;
        }
     } 
   } 
   free(CCb); free(X); 
}


static void splice_acceptor_uni(uchar *xseq, uchar *yseq, int M, int N,
                                int *aI_score, int **max_Ib, int **end_Ii)
{
   int *CCb, *X, *mI, *eI;
   int tmp, i, j, ss, ssx, cx, c;
   uchar *t, *s;

   CCb = (int *)ckalloc((M+1)*sizeof(int));
   X = (int *)ckalloc((M+1)*sizeof(int));
   mI = *max_Ib = (int *)ckalloc((2*N+2)*sizeof(int));
   eI = *end_Ii = mI+(N+1);

   t = yseq+N-1; CCb[M] = X[M] = 0;
   for (j=M-1; j>=0; j--) { CCb[j] = M-j; X[j] = 0; }
  
   mI[N] = -999999;        
   for (j=M; j>=0; j--)  
      if ((100*aI_score[j])>mI[N]) { mI[N] = 100*aI_score[j]; eI[N] = j+1; }
 
   for (i=N-1; i>=0; i--, t--) {
     s = xseq+M-1;
     ss = CCb[M]; ssx = X[M];
     c = ++CCb[M]; cx = X[M];
     for (j=M-1; j>=0; j--, s--) {
          tmp=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
          if (tmp==c+1) ;
          else if (tmp==CCb[j]+1) cx = X[j];
          else cx = ssx + (*t==*s);
 
          c = tmp; ss = CCb[j]; CCb[j] = c; ssx = X[j]; X[j] = cx;
     }
       
     /* compute max_Ib */           
     mI[i] = -999999; 
     for (j=M; j>=0; j--) {
        assert(CCb[j]+X[j]!=0);
        tmp = (int)(stepct(M-j)*X[j]/(double)(CCb[j]+X[j])*100)+100*aI_score[j];
        if (tmp>mI[i]) {
            mI[i] = tmp; eI[i] = j+1;
        }
     }
   }
   free(CCb); free(X);
}


static int stepct(int n)
{
  if (n<0) fatal("splice.c: Negative value in stepct().");
  if (n<=4) return 9;     
  if (n<=8) return 10;
  if (n<=12) return 12;
  return 12;
}

