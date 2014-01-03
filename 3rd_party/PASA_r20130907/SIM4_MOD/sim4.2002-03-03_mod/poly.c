#include "psublast.h"
#include "sim4.h"
#include "sim4b1.h"
#include "align.h" 
#include "poly.h"


#ifndef __lint
/*@unused@*/
static const char rcsid[] = "$Id: poly.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif

static void remove_polyT_front(struct edit_script_list **,Exon *,uchar *,uchar*,int *);
static void remove_polyA_back(struct edit_script_list **,Exon *,uchar *,uchar*,int,int *);
static void trim_polyT_align(struct edit_script_list **,Exon **,const int,int *,uchar *,uchar *);
static void trim_polyA_align(struct edit_script_list **,Exon *,Exon **,const int,int *,uchar *,uchar *);

void get_polyAT(uchar *seq, int len, int *pT, int *pA, int flag)
{
    register int i, sum10, sum20;
    register uchar *s, *t, *v;
    int last10;

    static char encodingA[128];
    static char encodingT[128];

    const int MAX10 = 2;
    const int MAX20 = 5;


    if (flag!=T_ONLY) {
       memset(encodingA, (char)1, 128);
       encodingA['A'] = encodingA['X'] = encodingA['N'] = 0;

       for (i=0, s=seq+len, sum10=0, last10=len+1; i<10 && s>seq && sum10<=MAX20; i++)  {
            sum10 += encodingA[*(--s)];
/*          if (!encodingA[*s] && sum10<=MAX10) last10 = s-seq+1; */
       }

       t = v = seq+len;
       sum20 = sum10;
       for ( ; s>=seq && (sum10<=MAX10 || sum20<=MAX20); ) {
            if (!encodingA[*s] && sum10<=MAX10 && (seq+len>=s+20 || sum20<MAX10))
                last10 = s-seq+1;
            if (--s>seq) {
                sum10 += encodingA[*s] - encodingA[*(--t)];
                sum20 += encodingA[*s] -(((seq+len)-s>20) ? encodingA[*(--v)] : 0);
            }
       }

       if (last10>len-10) *pA = len+1;
       else {
         s = seq+last10+8;
         while (*s && !encodingA[*s]) s--;
         if ((s-seq+1)-last10+1<=5)
             *pA = s-seq+2;
         else
             *pA = last10;
       }
    } else *pA = len+1;
    *pA = len-(*pA)+1;

    if (flag!=A_ONLY) {

       memset(encodingT, (char)1, 128);
       encodingT['T'] = encodingT['X'] = encodingT['N'] = 0;

       for (i=0, s=seq-1, sum10=0, last10=0; i<10 && i<len-1 && sum10<=MAX20; i++) {
            sum10 += encodingT[*(++s)];
/*          if (!encodingT[*s] && sum10<=MAX10) last10 = s-seq+1; */
       }

       t = v = seq-1;
       sum20 = sum10;
       for ( ; s<seq+len && (sum10<=MAX10 || sum20<=MAX20); ) {
            if (!encodingT[*s] && sum10<=MAX10 && (s-seq>=19 || sum20<MAX10))
                last10 = s-seq+1;
            if (++s<seq+len) {
                sum10 += encodingT[*s] - encodingT[*(++t)];
                sum20 += encodingT[*s] - ((s-seq>=20) ? encodingT[*(++v)] : 0);
            }
       }

       if (last10<=10) *pT = 0;
       else {
         s = seq+last10-10;
         while (*s && !encodingT[*s]) s++;
         if (last10-(s-seq)+1<=5)
             *pT = s-seq;
         else
             *pT = last10;
       }
    } else *pT = 0;
}


void remove_poly(struct edit_script_list **Script, Exon *Exons, uchar *s1, uchar *s2, int len2, int *pT, int *pA)
{
     remove_polyT_front(Script, Exons, s1, s2, pT); 
     remove_polyA_back(Script, Exons, s1, s2, len2, pA);
     *pA = len2-(*pA)+1;
/*   printf("pT: %d pA: %d\n", *pT, *pA); */

     return;
}

static void remove_polyA_back(struct edit_script_list **Sptr, Exon *Exons, uchar *s1, uchar *s2, int len2, int *lastA) {
     Exon *t, *exons_tail, *prev; /* start from Lblock */
     uchar *b, *end;
     int numA, pA, dummy, trim_p, reverse_script=0;

     *lastA = len2+1;  pA = 0;
     if (!Exons || ! Exons->next_exon || ! Exons->next_exon->to1) return;

     if ((*Sptr)->next_script &&
         (*Sptr)->offset1<(*Sptr)->next_script->offset1) {
         reverse_script = 1;
         script_flip_list(Sptr);
     }

     
     exons_tail = Exons->next_exon; prev = Exons;
     for ( ; exons_tail->next_exon && exons_tail->next_exon->to1; 
             prev=exons_tail, exons_tail=exons_tail->next_exon);

     trim_p = TRUE;
     while ((t=exons_tail)!=NULL && t->to1 && trim_p) {
        /* compute the 'A' contents of the exon */
        b = s2 + t->to2-1; end = s2+t->from2-1; numA = 0;
        while (b>=end && numA+(b-s2)>=.60*t->length) { 
               numA += (*b--=='A'); 
        }
        if (numA>=.60*t->length) {
            /* remove the entire exon */
            trim_polyA_align(Sptr,Exons,&exons_tail,t->from2,lastA,s1,s2);
/*          assert(*lastA==t->from2);  t was removed */
        } else {
            get_polyAT(s2+(*Sptr)->offset2-1,(*Sptr)->len2,&dummy,&pA,A_ONLY);
            if (pA) {
               int ct_pA;
               /* first position to be removed */ 
               ct_pA = t->to2-pA+1; 
               ct_pA = (ct_pA-t->from2>=MIN_EXON) ? ct_pA : t->from2;
               /* note: pA is the last (innermost) position in the tail */
               trim_polyA_align(Sptr,Exons,&exons_tail,ct_pA,lastA,s1,s2);
            }
            if (t==exons_tail) trim_p = FALSE;
        }
     }

     if (reverse_script) script_flip_list(Sptr);
}

static void trim_polyA_align(struct edit_script_list **Sptr, Exon *lblock, Exon **exons, const int bc, int *pA, uchar *s1,uchar *s2) 
{
    edit_script_list *head = *Sptr;
    edit_script *tp;
    int tmpi = 0, num, idents = 0;
    uchar *a, *b;
    Exon *prev;

    int i, j;  /* i index in the cDNA */

    if (bc>head->offset2+head->len2-1) {
        *pA = bc;
        return;
    }

    if (bc==head->offset2) {
        /* cDNA gap: remove the entire script; this should be properly sorted */
        *Sptr = head->next_script;
        Free_script(head->script);
        free(head);
        while ((*exons)->from2>=bc) {
           prev = find_previous(lblock,*exons);
           prev->next_exon = (*exons)->next_exon;
           free(*exons); *exons = prev;
        }
        *pA = bc;
        return;
    }

    Flip_script(&(head->script));
    i = head->offset2 + head->len2 -1;
    j = head->offset1 + head->len1 -1;
    tp = head->script;

    while (i>=bc && tp) {
       num = tp->num;
       switch (tp->op_type) {
          case INSERT:
                   if (i>=bc && bc>i-num+1) {
                       tmpi += i-bc+1; tp->num -= i-bc+1; i = bc-1;
                   } else {
                       i -= num; tmpi += num; head->script = tp->next;
                       free(tp); tp = head->script;
                   }
                   break;
          case DELETE:
                   j -= num; tmpi += num; head->script = tp->next;
                   free(tp); tp = head->script;
                   break;
          case SUBSTITUTE:
                   if (i>=bc && bc>i-num+1) {
                       a = s2+i-1; b = s1+j-1;
                       while (a>=s2+bc-1) {
                          if (*a--!=*b--) tmpi++; else idents++;
                       }
                       j -= i-bc+1; tp->num -= i-bc+1; i = bc-1;
                   } else {
                       /* at most 1 nt remaining */
                       a = s2+i-1; b = s1+j-1;
                       while (a>=s2+i-num) {
                          if (*a--!=*b--) tmpi++; else idents++;
                       }

                       i -= num; j -= num;
                       head->script = tp->next;
                       free(tp); tp = head->script;
                   }
                   break;
          default: fatalf("Unrecognized opcode %d.\n",tp->op_type);
       }
       /* indel walk */
    }
    assert(i==bc-1);

    while (tp->op_type!=SUBSTITUTE && j+1>=(*exons)->from1) {
       if (tp->op_type==INSERT) {
           i -= tp->num; tmpi += tp->num;
       } else if (j<(*exons)->from1 && i<(*exons)->from2) {
           j -= tp->num;
       } else {
           j -= tp->num; tmpi += tp->num;
       }
       head->script = tp->next;
       free(tp); tp = head->script;
    }

    if (head->script==NULL) {
        *Sptr = head->next_script;
        free(head);
    } else {
        head->len1 = j-head->offset1+1;
        head->len2 = i-head->offset2+1;
        head->score -= tmpi;
        Flip_script(&(head->script));
    }

    if ((*exons)->from2>i) {
        prev = find_previous(lblock,*exons);
        prev->next_exon = (*exons)->next_exon;
        free(*exons); *exons = prev;
    } else {
        double tmp_matches;
        (*exons)->to2 = i;
        (*exons)->to1 = j;
        (*exons)->length = (*exons)->to2-(*exons)->from2+1;
        tmp_matches = (*exons)->nmatches - idents;
        (*exons)->alen -= tmpi+idents;
        (*exons)->match = (int)(100*tmp_matches/(*exons)->alen);
    }
    *pA = i+1;

    return;
}


static void remove_polyT_front(struct edit_script_list **Sptr, Exon *Exons, uchar *s1, uchar *s2, int *lastT)
{
     Exon *t, *exons_head; /* start from Lblock */
     uchar *b, *end;
     int numT, dummy, trim_p, reverse_script=0, pT;

     *lastT = pT = 0;
     if (!Exons || !Exons->next_exon || !Exons->next_exon->to1) return;
 
     if ((*Sptr)->next_script && 
         (*Sptr)->offset1>(*Sptr)->next_script->offset1) {
         script_flip_list(Sptr);
         reverse_script = 1;
     }

     exons_head = Exons->next_exon; trim_p = TRUE;
     while ((t=exons_head)!=NULL && t->to1 && trim_p) {
        /* compute the 'T' contents of the exon */
        b = s2 + t->from2-1; end = s2+t->to2; numT = 0;
        while (b<end && (numT+t->to2-(b-s2-t->from2+1)>=.60*t->length)) {
                numT += (*b++=='T');
        }
        if (numT>=.60*t->length) {
            /* remove the entire exon */
            trim_polyT_align(Sptr,&exons_head,t->to2,lastT,s1,s2);
/*          assert(*lastT==t->to2);  t was removed */
        } else {
            get_polyAT(s2+(*Sptr)->offset2-1,(*Sptr)->len2,&pT,&dummy,T_ONLY);
            if (pT) {
                int ct_pT;
                ct_pT = pT + (*Sptr)->offset2-1;
                ct_pT = (t->to2-ct_pT>=MIN_EXON) ? ct_pT : t->to2;
                trim_polyT_align(Sptr,&exons_head,ct_pT,lastT,s1,s2);
            }
            if (t==exons_head) trim_p = FALSE;
        }
     }
     Exons->next_exon = exons_head;
     if (reverse_script) script_flip_list(Sptr);
}

/* s2 is the cdna */
static void trim_polyT_align(struct edit_script_list **Sptr, Exon **exons, const int ec, int *pT, uchar *s1, uchar *s2)
{
    edit_script_list *head = *Sptr;
    edit_script *tp;
    int tmpi = 0, num, idents = 0;
    uchar *a, *b;
    Exon *t;

    int i, j;  /* i index in the cDNA */

    if (ec<head->offset2) {
        *pT = ec; 
        return;
    }

    if (ec==head->offset2+head->len2-1) {
        /* cDNA gap: remove the entire script */
        *Sptr = head->next_script;
        Free_script(head->script);
        free(head);
        while ((*exons)->from2<ec) {
           t = *exons; *exons = t->next_exon; free(t);
        }
        *pT = ec;
        return;
    }
 
    i = head->offset2;
    j = head->offset1;
    tp = head->script;

    while (i<=ec && tp) {
       num = tp->num;
       switch (tp->op_type) {
          case INSERT:
                   if (i<=ec && ec<i+num-1) {
                       tmpi += ec-i+1; tp->num -= ec-i+1; i = ec+1; 
                   } else {
                       i += num; tmpi += num; head->script = tp->next; 
                       free(tp); tp = head->script; 
                   }
                   break;
          case DELETE:     
                   j += num; tmpi += num; head->script = tp->next;
                   free(tp); tp = head->script; 
                   break;
          case SUBSTITUTE:
                   if (i<=ec && ec<i+num-1) {
                       a = s2+i-1; b = s1+j-1;
                       while (a<s2+ec) {
                          if (*a++!=*b++) tmpi++; else idents++;
                       }
                       j += ec-i+1; tp->num -= ec-i+1; i = ec+1;
                   } else {
                       /* at most 1 nt remaining */
                       a = s2+i-1; b = s1+j-1;
                       while (a<s2+i+tp->num-1) {
                          if (*a++!=*b++) tmpi++; else idents++;
                       }

                       i +=num; j += num;
                       head->script = tp->next;
                       free(tp); tp = head->script; 
                   }
                   break;
          default: fatalf("Unrecognized opcode %d.\n",tp->op_type);
       }
       /* indel walk */
    }
    assert(i==ec+1);
     
    while (tp->op_type!=SUBSTITUTE && j-1<=(*exons)->to1) {
       if (tp->op_type==INSERT) {
           i += tp->num; tmpi += tp->num;
       } else if (j>=(*exons)->to1 && i>=(*exons)->to2) {
           j += tp->num;
       } else {
           j += tp->num; tmpi += tp->num;
       }
       head->script = tp->next;
       free(tp); tp = head->script;
    }
    
    if (head->script==NULL) {
        *Sptr = head->next_script;
        free(head);
    } else {
        head->len1 -= j-head->offset1;
        head->len2 -= i-head->offset2;
        head->offset2 = i;
        head->offset1 = j;          
        head->score -= tmpi;
    }

    if ((*exons)->to2<i) {
        t = *exons; *exons = t->next_exon; free(t);
    } else {
        double tmp_matches;
        (*exons)->from2 = i; 
        (*exons)->from1 = j;
        (*exons)->length = (*exons)->to2-(*exons)->from2+1;
        tmp_matches = (*exons)->nmatches - idents;
        (*exons)->alen -= tmpi+idents;
        (*exons)->match = (int)(100*tmp_matches/(*exons)->alen);
    }
    *pT = i-1;
    return;
}
