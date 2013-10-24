/*************************************************************************
*                                                                        *
*   Module: BuildDonors                                                  *
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

/*  $Id: BuildDonors.c,v 1.7 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=2 */
extern int TRANS[];

/* Maximum allowed number of generic sites */
extern long NUMSITES;
extern int UTR;



long  BuildDonors(char* s,
		  short class,
		  char* type,
		  char* subtype,
		  profile* p,
		  site* st, 
		  long l1, 
		  long l2,
		  long ns,
		  long nsites,
		  int Strand,
		  packExternalInformation* external
		  ) 
{
  int i,j;
  float score;
  long is;
  long left,right;
  int index;

  /* 1. Searching sites between beginning of the sequence and p->offset */
  if (!l1)
    {
      for (is = 0; is < p->offset && (ns<nsites); is++)
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
		  if(UTR){
		    score = score - PeakEdgeScore(is + p->order,Strand,external,l1,l2,6);
		  }
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = is + p->order;
			  st[ns].Score= score;
			  st[ns].class= class;
			  strcpy(st[ns].type,type);
			  strcpy(st[ns].subtype,subtype);
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
      while (*(s+p->dimension-1) && (is < right- left + 1) && (ns<nsites))
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
		  if(UTR){
		    score = score - PeakEdgeScore(left + is + p->offset,Strand,external,l1,l2,6);
		  }
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class= class;
			  strcpy(st[ns].type,type);
			  strcpy(st[ns].subtype,subtype);
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
      while (*(s+p->dimension-1) && (is < right- left + 1) && (ns<nsites))
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
		  if(UTR){
		    score = score - PeakEdgeScore(left + is + p->offset,Strand,external,l1,l2,6);
		  }
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class= class;
			  strcpy(st[ns].type,type);
			  strcpy(st[ns].subtype,subtype);
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
      while (*(s+p->dimension-1) && (is < right- left + 1) && (ns<nsites))
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
		  if(UTR){
		    score = score - PeakEdgeScore(left + is + p->offset,Strand,external,l1,l2,6);
		  }
		  if (score >= p->cutoff) 
			{
			  st[ns].Position = left + is + p->offset;
			  st[ns].Score=score;
			  st[ns].class= class;
			  strcpy(st[ns].type,type);
			  strcpy(st[ns].subtype,subtype);
			  ns++;
			}			
		  is++;
		  s++;
		}
    }
  
  if (ns >= nsites)
    printError("Too many predicted sites: decrease RSITES parameter");

  return(ns);
}

 
  
