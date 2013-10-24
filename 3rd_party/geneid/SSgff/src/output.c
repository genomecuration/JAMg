/*************************************************************************
*                                                                        *
*   Module: output                                                       *
*                                                                        *
*   Formatted output of SSgff                            )               *
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

extern int VRB;

/* Printing error messages */
void printError(char *s)
{
  fprintf(stderr,"Error: %s\n",s);
  exit(1);
}

/* Printing messages (information) */
void printMess(char* s)
{
  if (VRB) 
    fprintf (stderr, "> %s\n",s);
}


void printFasta (long pos1, long pos2, char strand, char *seq)
{
  long j,i=1; 
  /* For reverse strand */
  if (strand == '-')
	{
	  for(j = pos2 - 1; j >= pos1-1; j--,i++)
		{
		  printf("%c",complement(seq[j]));
		  if (!(i % 60))
			printf("\n");
		}
				 
	}
  else
	/* Forward  strand */
	{
	  for(j = pos1 - 1; j <= pos2 - 1; j++,i++)
		{
		  printf("%c", seq[j]);
		  if (!(i % 60))
			printf("\n");
		}
	  
	}
  if ((i-1) % 60)
	printf("\n"); 

}


void printHeadSite(char *Group, int numexon, char *nameseq,  char *Type)
{

  printf(">%s.%d:%s %s\n", Group, numexon+1, nameseq, Type);

}

void printHeadExon(char *Group, int numexon, char *nameseq)
{

  printf(">%s.%d:%s exon %d \n", Group, numexon+1, nameseq, numexon+1);

}

void printFastaGene(char *Group, char *nameseq, char *AuxGene)
{
  int j;

  printf(">%s:%s CDS\n", Group, nameseq);
  
  for(j=0; j < strlen(AuxGene); j++)
    {
	  printf("%c",AuxGene[j]);
      if (!((j+1)%60))
		printf("\n");
    }

  if (j %60)
	printf("\n"); 
}

void printHeadTranscript(char *Group, char *nameseq)
{

  printf(">%s:%s Primary Transcript \n", Group, nameseq);

}

void printHeadIntron(char *Group, int numexon, char *nameseq)
{

  printf(">%s.i%d:%s Intron between exon %d and %d \n", Group,  numexon, nameseq, numexon, numexon+1);
  
}





