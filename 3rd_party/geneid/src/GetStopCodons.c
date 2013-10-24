/*************************************************************************
*                                                                        *
*   Module: GetStopCodons                                                *
*                                                                        *
*   Stop codons prediction by using Position Weighted arrays             *
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

/*  $Id: GetStopCodons.c,v 1.10 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=2 */
extern int TRANS[];

/* Maximum allowed number of generic sites */
extern long NUMSITES;

long GetStopCodons(char* s,
                   profile* p,
                   site* sc, 
                   long l1, 
                   long l2) 
{
  long ns,is;
  float score;
  int i,j;
  
  /* Strings defining Stop codons */
  static char *stop[] =
  {
    "TAA",
    "TAG",
    "TGA"
  };

  char codon[LENGTHCODON+1];
  long left,right;
  int index;
  
  /* Final number of predicted stops */
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
		  
		  if (score >= p->cutoff) 
			{
			  sc[ns].Position = is + p->order;
			  sc[ns].Score=score;
			  sc[ns].class=U2;
			  ns++;
			}
		}
    }
  
  /* 2. Normal processing */
  /* left and right are the true boundaries of prediction */
  left  = MAX(0+p->order,l1 - p->offset);
  right = l2 - p->offset;
  s += left;
  is = 0;     
  
  /* Case A: Using Markov chain with order 0: PWM */
  if (p->order==0)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
		{
		  /* The candidate region must contain the STOP in the right position */
		  strncpy(codon,s+p->offset+1,LENGTHCODON);
		  codon[LENGTHCODON]='\0';
		  
		  if (!strcmp(codon,stop[0]) ||
			  !strcmp(codon,stop[1]) ||
			  !strcmp(codon,stop[2])) 
			{
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
			  
			  if (score >= p->cutoff) 
				{
				  /* Position given is last coding position before Stop */ 
				  sc[ns].Position=left + is + p->offset;
				  sc[ns].Score=score;
				  sc[ns].class=U2;
				  ns++;
				}
			} 
		  is++;
		  s++;
		}
    }   
  /* Case B: Using Markov chain with order 1: dinucleotides */
  else if (p->order==1)
    {
      /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
		{
		  /* The candidate region must contain the STOP in the right position */
		  strncpy(codon,s+p->offset+1,LENGTHCODON);
		  codon[LENGTHCODON]='\0';
		  
		  if (!strcmp(codon,stop[0]) ||
			  !strcmp(codon,stop[1]) ||
			  !strcmp(codon,stop[2])) 
			{
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
			  
			  if (score >= p->cutoff) 
				{
				  /* Position given is last coding position before Stop */ 
				  sc[ns].Position=left + is + p->offset;
				  sc[ns].Score=score;
				  sc[ns].class=U2;
				  ns++;
				}
			} 
		  is++;
		  s++;
		}     
    }
  /* case C: Using Markov chain with order > 1 */
  else
    {
	  /* discovering splice sites with current profile */
      while (*(s+p->dimension) && (is < right - left + 1) && (ns<NUMSITES))  
		{
		  /* The candidate region must contain the STOP in the right position */
		  strncpy(codon,s+p->offset+1,LENGTHCODON);
		  codon[LENGTHCODON]='\0';
		  
		  if (!strcmp(codon,stop[0]) ||
			  !strcmp(codon,stop[1]) ||
			  !strcmp(codon,stop[2])) 
			{
			  score=0.0;
			  for (i=0;i<p->dimension;i++) 
				{
				  /* i is the position inside the region */
				  index = OligoToInt(s + i - p->order , p->order+1,5); 
				  if (index >= p->dimensionTrans)
					score = score + -INFI;
				  else
					score = score + p->transitionValues[i][index];
				}
			  
			  if (score >= p->cutoff) 
				{
				  /* Position given is last coding position before Stop */ 
				  sc[ns].Position=left + is + p->offset;
				  sc[ns].Score=score;
				  sc[ns].class=U2;
				  ns++;
				}
			} 
		  is++;
		  s++;
		} 
    }
  
  /* 3. Remaining stops until the end of sequence. Set score to 0 for those */
  if (!(s+p->dimension))
    {
      s=(s-is);
      is+=p->offset;
      
      while (*(s+is) && (ns<NUMSITES))  
		{
		  /* The candidate region must contain the STOP in the right position */
		  strncpy(codon,(s+is),LENGTHCODON);
		  codon[LENGTHCODON]='\0';
		  
		  if (!strcmp(codon,stop[0]) || !strcmp(codon,stop[1]) || !strcmp(codon,stop[2]))
			{
			  /* Position given is last coding position before Stop */
			  sc[ns].Position=is-1;  
			  sc[ns].Score=0;
			  sc[ns].class=U2;
			  ns++;
			} 
		  is++;
		}
    }
  
  if (ns >= NUMSITES)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);
}
