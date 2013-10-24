/*************************************************************************
*                                                                        *
*   Module: GetSitesWithProfile                                          *
*                                                                        *
*   Signal prediction by using a Position Weighted Array                 *
*                                                                        *
*   This file is part of the geneid 1.4 distribution                     *
*                                                                        *
*     Copyright (C) 2006 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           *
*                          Tyler   ALIOTO                                * 
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           * 
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

/*  $Id: GetSitesWithProfile.c,v 1.11 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=2 */
extern int TRANS[];

/* Maximum allowed number of generic sites */
extern long NUMSITES;

long  GetSitesWithProfile(char* s,
                          profile* p,
                          site* st, 
                          long l1, 
                          long l2) 
{
  int i,j;
  float score;
  long ns,is;
  long left,right;
  int index;
  
  /* Final number of predicted signals (that type) */
  ns = 0;
  
  /* 1. Searching sites between beginning of the sequence and p->offset */
  if (!l1)
    {
      for (is = 0; is < p->offset && (ns<NUMSITES); is++)
		{
		  score=0.0;
		  /* Applying part of the profile */
		  for (i=p->offset-is, j=0; i < p->dimension; i++,j++) 
			{
			  /* i is the position inside the region */
			  index = OligoToInt(s+j, p->order+1,5);
			  
			  if (index >= p->dimensionTrans)
				score = score + -INFI;
			  else
				score = score + p->transitionValues[i][index];
			}
        
		  score = p->afactor + (p->bfactor * score);
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = is + p->order;
			  st[ns].Score=score;
			  st[ns].class=U2;
/* 			  strcpy(st[ns].subtype,sU2); */
			  ns++;
			}
		}
    }
  
  /* 2. Normal processing: predicting using the whole profile */
  /* left and right are the true boundaries of prediction */
  left  = MAX(0+p->order, l1 - p->offset);
  right = l2 - p->offset;
  s += left;
  is = 0;     
  /* Case A: Using Markov chain with order 0: PWM */
  if (p->order == 0)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right- left + 1) && (ns<NUMSITES))
		{ 
		  /* is = 0..right */
		  score=0.0;
		  for (i=0;i<p->dimension;i++)
			{
			  /* i is the position inside the region */
			  index = TRANS[(int)(*(s + i))];
			  if (index >= p->dimensionTrans)
				score = score + -INFI;
			  else
				score = score + p->transitionValues[i][index];
			}
		  
		  score = p->afactor + (p->bfactor * score);
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class=U2;
/* 			  strcpy(st[ns].subtype,sU2); */
			  ns++;
			}
		  is++;
		  s++;
		}
    }
  /* case B: Using Markov chain with order 1: dinucleotides */
  else if (p->order == 1)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right- left + 1) && (ns<NUMSITES))
		{ 
		  /* is = 0..right */
		  score=0.0;
		  for (i=0;i<p->dimension;i++)
			{
			  /* i is the position inside the region */
			  index = 5*TRANS[(int)(*(s + i -1))] + TRANS[(int)(*(s + i))];
			  if (index >= p->dimensionTrans)
				score = score + -INFI;
			  else
				score = score + p->transitionValues[i][index];
			}
		  
		  score = p->afactor + (p->bfactor * score);
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class=U2;
/* 			  strcpy(st[ns].subtype,sU2); */
			  ns++;
			}
		  is++;
		  s++;
		}
    }
  /* case C: Using Markov chain with order > 1 */
  else
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right- left + 1) && (ns<NUMSITES))
		{ 
		  /* is = 0..right */
		  score=0.0;
		  for (i=0;i<p->dimension;i++)
			{
			  /* i is the position inside the region */
			  /* 5 is used because there are A,C,G,T and N */
			  index = OligoToInt(s + i - p->order, p->order+1,5);

			  if (index >= p->dimensionTrans)
				score = score + -INFI;
			  else
				score = score + p->transitionValues[i][index];
			}
		  
		  score = p->afactor + (p->bfactor * score);
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class=U2;
/* 			  strcpy(st[ns].subtype,sU2); */
			  ns++;
			}
		  is++;
		  s++;
		}
    }
  
  if (ns >= NUMSITES)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);
}

 
  
