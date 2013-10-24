/*************************************************************************
*                                                                        *
*   Module: scan                                                         *
*                                                                        *
*   Formatted output of geneid (GFF, default and extended)               *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Genis PARRA FARRE                             *
*                          Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
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

#include "SSgff.h"

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
		// preliminary change before changing readseq 
      case 'a':
        return('t');
      case 'c':
        return('g');
      case 'g':
        return('c');
      case 't':
        return('a');

      default:
        return('N');
      }
}

/* Reverse a split of the original sequence */
void ReverseSubSequence (long p1, long p2, char* s, char* r)
{
  long i,j;
  long l;

  l = p2 - p1 + 1;
  for (i = p1, j = p1 + (l-1); i <= j; i++, j--)
    {    
      /* Reverse */
      r[i-p1]=complement(s[j]);
      r[j-p1]=complement(s[i]);
    }
}

void CheckBoundaries (long *pos1, long *pos2, long length)
{
  if (*pos1 < 1)
	{	
	  *pos1 = 1;
	  printMess("WARNING!: coordinate out of the range");
	}
  if (*pos2 > length)
	{
	  *pos2 = length;
	  printMess("WARNING!: coordinate out of the range");
	}
}


void ScanSequence (long pos1, long pos2, char *Sequence, char strand, char *saux)
{

  /* Extrancting the subsequence */
  if (strand == '-')
	/* Reverse Strand */
    ReverseSubSequence(pos1 - 1, pos2 - 1, Sequence, saux);
  else 
	{ 
	  /* Forward Strand */  
	  strncpy(saux, Sequence + pos1 - 1, pos2 - pos1 +1);
	}

  /* including end string */
  saux[pos2 - pos1 + 1] = '\0';

}


