/*************************************************************************
*                                                                        *
*   Module: FetchSequence                                                *
*                                                                        *
*   Reverse and complement DNA sequences (upper and lower case)          *
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

/*  $Id: FetchSequence.c,v 1.5 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Complement a nucleotide */
int complement(int c)
{
  switch (c) {
  case 'A':
	return('T');
  case 'C':
	return('G');
  case 'G':
	return('C');
  case 'T':
	return('A');
  default:
	return('N');
  }
}

/* Translate into upper letters and reverse/complement a DNA sequence */
/* s is the input and r is the output (reverse) */
long FetchSequence(char *s, char* r)
{
  long i,j;
  long l;
  
  /* Length of input sequence */
  l = strlen(s);
  
  /* Scan input sequence, producing the reversed and complemented output */
  /* Run from both ends of input sequence in a parallel way */
  for (i = 0, j = l-1; i <= j; i++, j--)
    {    
      /* Upper case */
      if (s[i] > 96)
		s[i] = s[i] - 32;
      if (s[j] > 96)
		s[j] = s[j] - 32;
	  
      /* Reverse(i,j) and complement */
      r[i]=complement(s[j]);
      r[j]=complement(s[i]);
    }
  
  r[l] = '\0';
  
  return (l);
}

/* Reverse a split in the original sequence (p1,p2) */
void ReverseSubSequence(long p1, long p2, char* s, char* r)
{
  long i,j;
  long l;
  
  /* Length of the input subsequence (fragment) */
  l = p2 - p1 + 1;
  
  /* Scan input fragment, producing the reversed and complemented output */
  for (i = p1, j = p1 + (l-1); i <= j; i++, j--)
    {    
      /* Reverse */
      r[i-p1]=complement(s[j]);
      r[j-p1]=complement(s[i]);
    }
}


