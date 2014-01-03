#include <stdio.h>
#include <stdlib.h>

#include "psublast.h"

#include "sim4.h"
#include "sim4b1.h"
#include "Xtend1.h"
#include "align.h"

#ifndef __lint
/*@unused@*/
static const char rcsid[] =
"$Id: align.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

static int snake(int k, int x, int endx, int endy);
static int rsnake(int k, int x, int startx, int starty, int M);


void    align_path(int i1, int j1, int i2, int j2, int dist, edit_script **head, edit_script **tail)
{

        int     *last_d, *temp_d,       /* forward vectors */
                *rlast_d, *rtemp_d;     /* backward vectors */

        edit_script *head1, *tail1, *head2, *tail2;
        int midc, rmidc;
        int start, lower, upper;
        int rstart, rlower, rupper;
        int c, k, row;
        int mi, mj, tmp, ll, uu;
        char flag;

        *head = *tail = NULL;

        /* Boundary cases */
        if (i1 == i2) {
           if (j1 == j2) *head = NULL;
           else {
                head1 = (edit_script *) ckalloc(sizeof(edit_script));
                head1->op_type = INSERT;
                head1->num = j2-j1;
                head1->next = NULL;
                *head = *tail = head1;
           }
           return;
        }

        if (j1 == j2) {
                head1 = (edit_script *) ckalloc(sizeof(edit_script));
                head1->op_type = DELETE;
                head1->num = i2-i1;
                head1->next = NULL;
                *head = *tail = head1;
                return;
        }

        if (dist <= 1) {
           start = j1-i1; 
           if (j2-i2 == j1-i1) {
                head1 = (edit_script *) ckalloc(sizeof(edit_script));
                head1->op_type = SUBSTITUTE;
                head1->num = i2-i1;
                head1->next = NULL;
                *head = *tail = head1;
           } else if (j2-j1 == i2-i1+1) {

                tmp = snake(start,i1,i2,j2);
                if (tmp>i1) {
                    head1 = (edit_script *) ckalloc(sizeof(edit_script));
                    head1->op_type = SUBSTITUTE;
                    head1->num = tmp-i1;
                    *head = head1;
                }
                head2 = (edit_script *) ckalloc(sizeof(edit_script));
                head2->op_type = INSERT;
                head2->num = 1;

                if (*head) head1->next = head2; 
                else *head = head2; 
                *tail = head2;
                head2->next = NULL;

                if (i2-tmp) {
                    head1 = head2;
                    *tail = head2 = (edit_script *)ckalloc(sizeof(edit_script));
                    head2->op_type = SUBSTITUTE;
                    head2->num = i2-tmp;
                    head2->next = NULL;
                    head1->next = head2;
                } 
           } else if (j2-j1+1 == i2-i1) {

                tmp = snake(start,i1,i2,j2);
                if (tmp>i1) {
                    head1 = (edit_script *) ckalloc(sizeof(edit_script));
                    head1->op_type = SUBSTITUTE;
                    head1->num = tmp-i1;
                    *head = head1;
                }
                head2 = (edit_script *) ckalloc(sizeof(edit_script));
                head2->op_type = DELETE;
                head2->num = 1;

                if (*head) head1->next = head2;
                else *head = head2;
                *tail = head2;
                head2->next = NULL;

                if (i2>tmp+1) {
                    head1 = head2;
                    *tail = head2 = (edit_script *)ckalloc(sizeof(edit_script));
                    head2->op_type = SUBSTITUTE;
                    head2->num = i2-tmp-1;
                    head2->next = NULL;
                    head1->next = head2;
                }
           } else {
                (void)fprintf(stderr, 
                      "align.c: warning: something wrong when aligning."); 
           }                
           return;
        }
                
        /* Divide the problem at the middle cost */
        midc = dist/2;
        rmidc = dist - midc;
                
        /* Compute the boundary diagonals */
        start = j1 - i1;
        lower = max(j1-i2, start-midc);
        upper = min(j2-i1, start+midc);
        rstart = j2-i2; 
        rlower = max(j1-i2, rstart-rmidc);
        rupper = min(j2-i1, rstart+rmidc);

        /* Allocate space for forward vectors */
        last_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
        temp_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;

        for (k=lower; k<=upper; k++) last_d[k] = -1;
        last_d[start] = snake(start,i1,i2,j2);
                
        /* Forward computation */  
        for (c=1; c<=midc; ++c) {
             ll = max(lower,start-c); 
             uu = min(upper,start+c);
             for (k=ll; k<=uu; ++k) {
                  if (k == ll) {
                      /* DELETE : down from (k+1,c-1) */
                      row = last_d[k+1]+1;
                  } else if (k == uu) {
                      /* INSERT : right from (k-1,c-1) */
                      row = last_d[k-1];
                  } else if ((last_d[k]>=last_d[k+1]) &&
                             (last_d[k]+1>=last_d[k-1])) {
                      /* SUBSTITUTE */
                      row = last_d[k]+1;
                  } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                             (last_d[k+1]>=last_d[k])) {
                      /* DELETE */
                      row = last_d[k+1]+1;
                  } else {
                      /* INSERT */
                      row = last_d[k-1];
                  }

                  temp_d[k] = snake(k,row,i2,j2);
             } 
             for (k=ll; k<=uu; ++k)
                  last_d[k] = temp_d[k];
       }

        /* Allocate space for backward vectors */
        rlast_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;
        rtemp_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;

        for (k=rlower; k<=rupper; k++) rlast_d[k] = i2+1;
        rlast_d[rstart] = rsnake(rstart,i2,i1,j1,M);

        /* Backward computation */
        for (c=1; c<=rmidc; ++c) {
             ll = max(rlower,rstart-c);
             uu = min(rupper,rstart+c);
             for (k=ll; k<=uu; ++k) {
                  if (k == ll) {
                      /* INSERT : left from (k+1,c-1) */
                      row = rlast_d[k+1];
                  } else if (k == uu) {
                      /* DELETE : up from (k-1,c-1) */  
                      row = rlast_d[k-1]-1;
                  } else if ((rlast_d[k]-1<=rlast_d[k+1]) &&
                             (rlast_d[k]-1<=rlast_d[k-1]-1)) {
                      /* SUBSTITUTE */
                      row = rlast_d[k]-1;
                  } else if ((rlast_d[k-1]-1<=rlast_d[k+1]) &&
                             (rlast_d[k-1]-1<=rlast_d[k]-1)) {
                      /* DELETE */
                      row = rlast_d[k-1]-1;
                  } else {
                      /* INSERT */
                      row = rlast_d[k+1];
                  }
               
                  rtemp_d[k] = rsnake(k,row,i1,j1,M);
             }    
             for (k=ll; k<=uu; ++k)
                  rlast_d[k] = rtemp_d[k];
       }

       /* Find (mi, mj) such that the distance from (i1, j1) to (mi, mj) is
          midc and the distance from (mi, mj) to (i2, j2) is rmidc.
        */

       flag = FALSE;
       mi = i1; mj = j1;
       ll = max(lower,rlower); uu = min(upper,rupper);
       for (k=ll; k<=uu; ++k) {
            if (last_d[k]>=rlast_d[k]) {
                if (last_d[k]-i1>=i2-rlast_d[k]) {
                    mi = last_d[k]; mj = k+mi;
                } else {
                    mi = rlast_d[k]; mj = k+mi; 
                } 
                flag = TRUE;

                break;
            }         
       }              
       free(last_d+lower); free(rlast_d+rlower);
       free(temp_d+lower); free(rtemp_d+rlower);
                      
       if (flag) {    
                /* Find a path from (i1,j1) to (mi,mj) */
                align_path(i1,j1,mi,mj,midc,&head1,&tail1);
                      
                /* Find a path from (mi,mj) to (i2,j2) */
                align_path(mi,mj,i2,j2,rmidc,&head2,&tail2);
                      
                /* Join these two paths together */
                if (head1) tail1->next = head2;
                else head1 = head2;
        } else {  
                (void)fprintf(stderr, 
                      "align.c: warning: something wrong when dividing\n");
                head1 = NULL;
        }         
        *head = head1;
        if (head2) *tail = tail2;
        else *tail = tail1;
}


int align_get_dist(int i1, int j1, int i2, int j2, int limit)
{
        int *last_d, *temp_d;
        int goal_diag, ll, uu;
        int c, k, row;
        int start, lower, upper;
                
        /* Compute the boundary diagonals */
        start = j1 - i1;
        lower = max(j1-i2, start-limit);
        upper = min(j2-i1, start+limit);
        goal_diag = j2-i2;
       
        if (goal_diag > upper || goal_diag < lower) {
           (void)fprintf(stderr, "The two sequences are not really similar.\n");
           (void)fprintf(stderr, "Please try an exact aligning method.\n");
           exit(1);
        }       
                      
        /* Allocate space for forward vectors */ 
        last_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
        temp_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;

        /* Initialization */
        for (k=lower; k<=upper; ++k) last_d[k] = MININT;
        last_d[start] = snake(start, i1, i2, j2);
        
        if (last_d[goal_diag] >= i2) {
                /* Free working vectors */
                free(last_d+lower);
                free(temp_d+lower);
                return 0;
        }

        for (c=1; c<=limit; ++c) {
          ll = max(lower,start-c); uu = min(upper, start+c);
          for (k=ll; k<=uu; ++k) {
               if (k == ll)
                        row = last_d[k+1]+1;    /* DELETE */
               else if (k == uu)
                        row = last_d[k-1];      /* INSERT */
               else if ((last_d[k]>=last_d[k+1]) &&
                             (last_d[k]+1>=last_d[k-1]))
                        row = last_d[k]+1;      /*SUBSTITUTE */
               else if ((last_d[k+1]+1>=last_d[k-1]) &&
                             (last_d[k+1]>=last_d[k]))
                        row = last_d[k+1]+1;    /* DELETE */
               else
                        row = last_d[k-1];      /* INSERT */
                
               temp_d[k] = snake(k,row,i2,j2);
          }     
                      
          for (k=ll; k<=uu; ++k) last_d[k] = temp_d[k];

          if (last_d[goal_diag] >= i2) {
#ifdef STATS
                 (void)fprintf(stderr, "get_dist = %d\n",c);
#endif
        
                 /* Free working vectors */
                    free(last_d+lower);
                    free(temp_d+lower);
                    return c;
           }
        }

        /* Ran out of distance limit */
        return -1;
}

/* Condense_script - merge contiguous operations of the same type together */
void Condense_script(edit_script *head)
{
        edit_script *tp, *tp1;

        tp = head;
        while (tp != NULL) {
           while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
                 tp->num = tp->num + tp1->num;
                 tp->next = tp1->next;
                 free(tp1);
           }
           tp = tp->next;
        }
}

/* Condense_both_Ends  --  merge contiguous operations of the same type    */
/* together; return both new ends of the chain.                            */
void Condense_both_Ends (edit_script **head, edit_script **tail, edit_script **prev)
{
        edit_script *tp, *tp1;

        tp = *head; *prev = NULL;
        while (tp != NULL) {
           while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
                 tp->num = tp->num + tp1->num;
                 tp->next = tp1->next;
                 free(tp1);
           }
           if (tp->next) *prev = tp;
           else *tail = tp;
           tp = tp->next;
        }
}

/* Flip_script - reverse the script list */
void Flip_script(struct edit_script **script)
{
   struct edit_script *ep, *ahead, *behind;

   ahead = *script;
   ep = NULL;
   while (ahead!=NULL) {
          behind = ep;
          ep = ahead;
          ahead = ahead->next;
          ep->next = behind;
  }
  *script = ep;
}

/* reverse the edit script */
#ifdef AUXUTILS
void Reverse_script(edit_script *head)
{
        edit_script *tp;

        tp = head;
        while (tp != NULL) {
           if (tp->op_type == INSERT) tp->op_type = DELETE;
           else if (tp->op_type == DELETE) tp->op_type = INSERT;
           tp = tp->next;
        }
}
#endif


#ifdef AUXUTILS
void Print_script(edit_script *head, int M, int N)
{
        edit_script *tp;
        int i, j, k, count;

        i = j = 0;
        tp = head;

        while (tp != NULL) {
           if (tp->op_type == SUBSTITUTE) {
                k = 0;
                while (k < tp->num) {
                        count = 0;
                        while ((seq1[i] == seq2[j]) &&
                               (k<tp->num)) {
                                ++i;
                                ++j;
                                ++count;
                                ++k;
                        }
                        if (count > 0) (void)printf("copy %d\n", count);
                        if (k < tp->num) {
                                (void)printf("replace %c by %c\n", seq1[i], seq2[j]);
                                ++i;
                                ++j;
                                ++k;
                        }
                }
                /*
                if (tp->num > 1) (void)printf("%d substitutions\n", tp->num);
                else (void)printf("%d substitution\n", tp->num);
                */
           } else if (tp->op_type == INSERT) {
                if ((tp==head || tp->next==NULL) && (M <= N))
                        (void)printf("skip (second sequence) %d\n", tp->num);
                else {
                        /*
                        (void)printf("insert %d\n", tp->num);
                        */
                        (void)printf("insert ");
                        if (tp->num<=10)
                            for (k=j; k<j+tp->num; ++k) (void)printf("%c", seq2[k]);
                        else (void)printf(" %d ",tp->num);
                        (void)printf("\n");
                }
                j += tp->num;
           } else {     /* DEL */
                if ((tp==head || tp->next==NULL) && (M > N))
                        (void)printf("skip (first sequence) %d\n", tp->num);
                else {
                        /*
                        (void)printf("delete %d\n", tp->num);
                        */
                        (void)printf("delete ");
                        if (tp->num <= 10)
                            for (k=i; k<i+tp->num; ++k)
                                 (void)printf("%c", seq1[k]);
                        else
                            (void)printf("%d ",tp->num);
                        (void)printf("\n");
                }
                i += tp->num;
           }
           tp = tp->next;
        }
}
#endif

void S2A(edit_script *head, int *S, int flag)
{
        edit_script *tp;
        int *lastS, i;

        tp = head;
        lastS = S;
        while (tp != NULL) {
/*
        (void)printf("tp->op_type=%d, tp->num=%d\n",tp->op_type, tp->num);
*/
           if (tp->op_type == SUBSTITUTE) {
                for (i=0; i<tp->num; ++i) *lastS++ = 0;
           } else if (tp->op_type == INSERT) {
                *lastS++ = (!flag) ? tp->num : (0-tp->num);
           } else {     /* DELETE */
                *lastS++ = (!flag) ? (0 - tp->num) : tp->num;
           }
           tp = tp->next;
        }
        *(S-1) = lastS - S;
}

void align_reverse(int *S)
{
       int auxi, *begi, *endi;
        
       begi = S; endi = S + *(S-1);
       while (begi < endi) {
          auxi = *begi;
          *begi = *--endi;
          *endi = auxi;
          begi++;      
       }               
       return;         
}                
                         
/* Alignment display routine */

static uchar ALINE[51], BLINE[51], CLINE[51];

void IDISPLAY(uchar A[], uchar B[], int M, int N, int S[], int AP, int BP,
              int est_strand, Exon *exons)
{ 
  Exon *t0;
  register uchar *a, *b, *c, sign;
  register int    i,  j, op, index;
  int   lines, ap, bp, starti;
        
  if ((exons==NULL) || (!exons->to1 && (exons->next_exon==NULL)))
       fatal("align.c: Exon list cannot be empty.");
           
  /* find the starting exon for this alignment */
  t0 = exons;
  while (t0 && (((est_strand==2) && ((t0->from1!=AP) || (t0->from2!=BP))) ||
                ((est_strand==1) && ((t0->from1!=BP) || (t0->from2!=AP)))))
     t0 = t0->next_exon;
  if (!t0) fatal("align.c: Alignment fragment not found.");
          
  i = j = op = lines = index = 0;
  sign = '*'; ap = AP; bp = BP; a = ALINE; b = BLINE; c = CLINE;
  starti = (t0->next_exon && t0->next_exon->to1) ? (t0->to1+1):-1;

  while (i < M || j < N) {
    if (op == 0 && *S == 0) {
          op = *S++; *a = A[++i]; *b = B[++j];
          *c++ = (*a++ == *b++) ? '|' : ' ';
    } else {
        if (op == 0) { op = *S++; }
        if (op > 0) {
           if (est_strand==2) {
               *a++ = ' '; *b++ = B[++j]; *c++ = '-'; op--;
           } else {
               if (j+BP==starti) {
                   /* detected intron */
                   switch (t0->ori) {
                      case 'C': sign = '<'; break;
                      case 'G': sign = '>'; break;
                      case 'N': sign = '='; break;
                      default: fatal("align.c: Unrecognized intron type.");
                   } 
                   t0 = t0->next_exon;
                   starti=(t0->next_exon && t0->next_exon->to1)?(t0->to1+1):-1;
                   index = 1; *c++ = sign; *a++ = ' '; *b++ = B[++j]; op--;
               } else if (!index) {
                   *c++ = '-'; *a++ = ' '; *b++ = B[++j]; op--;
               } else { 
                   /* not the first deletion in the intron */
                   switch (index) {
                       case 0:
                       case 1:
                       case 2: *a++ = ' '; *b++ = B[++j];
                               *c++ = sign; op--; index++; break;
                       case 3:
                       case 4: *a++ = ' '; *b++ = '.'; *c++ = '.';
                                j++; op--; index++; break;
                       case 5: *a++ = ' '; *b++ = '.'; *c++ = '.';
                                j+= op-3; op = 3; index++; break;
                       case 6:
                       case 7: *a++ = ' '; *b++ = B[++j]; *c++ = sign;
                               op--; index++; break;
                       case 8: *a++ = ' '; *b++ = B[++j];
                               *c++ = sign; op--; index = 0; break;
                       }
               }   
           }   
        } else {   
           if (est_strand==1) {
               *a++ = A[++i]; *b++ = ' '; *c++ = '-'; op++;  
           } else {
               if (i+AP==starti) {
                   /* detected intron */
                   switch (t0->ori) { 
                      case 'C': sign = '<'; break;
                      case 'G': sign = '>'; break;
                      case 'N': sign = '='; break;
                      default: fatal("align.c: Unrecognized intron type.");
                   }   
                   t0 = t0->next_exon;
                   starti=(t0->next_exon && t0->next_exon->to1)?(t0->to1+1):-1;
                       
                   index = 1; *c++ = sign; *a++ = A[++i]; *b++ = ' '; op++;
               } else if (!index) { 
                   *c++ = '-'; *a++ = A[++i]; *b++ = ' '; op++;
               } else { 
                   /* not the first deletion in the intron */
                   switch (index) {
                       case 0:
                       case 1: 
                       case 2: *a++ = A[++i]; *b++ = ' '; *c++ = sign; op++;
                                index++; break;
                       case 3:
                       case 4: *a++ = '.'; *b++ = ' '; *c++ = '.';
                                i++; op++; index++; break;
                       case 5: *a++ = '.'; *b++ = ' '; *c++ = '.';
                                i+=(-op)-3; op=-3; index++; break;
                       case 6:
                       case 7: *a++ = A[++i]; *b++ = ' '; 
                               *c++ = sign; op++; index++; break;
                       case 8: *a++ = A[++i]; *b++ = ' ';
                               *c++ = sign; op++; index = 0; break;
                   }   
               }   
           }   
        }          
    }          
    if ((a >= ALINE+50) || ((i >= M) && (j >= N))) {
        *a = *b = *c = '\0';
        (void)printf("\n%7d ",50*lines++);
        for (b = ALINE+10; b <= a; b += 10)
            (void)printf("    .    :");
        if (b <= a+5)           
            (void)printf("    .");
        (void)printf("\n%7d %s\n        %s\n%7d %s\n",ap,ALINE,CLINE,bp,BLINE);
         ap = AP + i;           
         bp = BP + j;  
         a = ALINE;             
         b = BLINE;    
         c = CLINE;
    }
  }
}


void Free_script(edit_script *head)
{
        edit_script *tp, *tp1;

        tp = head;
        while (tp != NULL) {
           tp1 = tp->next;
           free(tp);
           tp = tp1;
        }
}


static int snake(int k, int x, int endx, int endy)
{
        int y;

        if (x<0) return x;
        y = x+k;
        while (x<endx && y<endy && seq1[x]==seq2[y]) {
                ++x; ++y;
        }
        return x;
}


static int rsnake(int k, int x, int startx, int starty, int M)
{
        int y;

        if (x>M) return x;
        if ((startx<0) || (starty<0))
            (void)printf("TROUBLE!!! startx:  %5d,  starty:  %5d\n",startx, starty);
        if ((x>M) || (x+k>N))
             (void)printf("TROUBLE!!! x:  %5d,  y:  %5d\n",x,x+k);
        
        y = x+k;
        while (x>startx && y>starty && seq1[x-1]==seq2[y-1]) {
                --x; --y;
        }
        return x;
}


