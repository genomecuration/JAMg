/*************************************************************************
*                                                                        *
*   Module: BuildORFs                                                    *
*                                                                        *
*   From start and stop codons, to build ORFs                            *
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

/*  $Id: BuildORFs.c,v 1.5 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (divided by RSINGL) */
extern long NUMEXONS;
extern long MAXBACKUPSITES;

long BuildORFs(site *Start, long nStarts, 
               site *Stop, long nStops,
               long cutPoint, char* Sequence,
               exonGFF *Exon) 
{
  int Frame;
  long i, j, js;
  
  /* Maximum allowed number of ORFs per fragment */     
  long HowMany;
  
  /* Final number of predicted ORFs */
  long nSingles;
  
  /* Main loop, for each Start codon searching the first Stop in frame */
  HowMany=(MAXBACKUPSITES)? (long)(NUMEXONS/RSINGL): NUMEXONS;
  for (i=0, j=0, nSingles=0;
	   (i < nStarts) && (j<nStops) && (nSingles<HowMany);
	   i++)
    {
      Frame = ((Start+i)->Position + 1) % 3;
	  
      /* Skip previous Stops to Start */
      while ( (j < nStops) && (((Stop+j)->Position+1) <= (Start+i)->Position+1))
		j++;
	  
      /* Save counter j for the next iteration */
      js=j;
	  
      /* Skip Stops not in frame with the current Start */
      while ((js < nStops) && (((Stop+js)->Position+1) % 3 != Frame))
		js++;
	  
      /* CutPoint: to preserve sorted exons between fragments */
      if (js < nStops && (Stop+js)->Position >= cutPoint)
		{
		  /* LENGTH rule about ORFs */
		  if ( ((Stop+js)->Position + LENGTHCODON - (Start+i)->Position + 1)
			   >= 
			   ORFLENGTH)
			{
			  (Exon + nSingles)->Acceptor = (Start+i);
			  (Exon + nSingles)->Donor = (Stop+js);
			  (Exon + nSingles)->Frame = 0;
			  (Exon + nSingles)->Remainder = 0;
			  strcpy((Exon + nSingles)->Type,"ORF");
			  strcpy((Exon + nSingles)->Group,NOGROUP);
			  (Exon + nSingles)->evidence = 0;

			  /* Store info about frame and remainder nucleotides to avoid building stops in frame */
			  ComputeStopInfo((Exon+nSingles),Sequence);
		 
			  nSingles++;
			}
		}
    }
  
  if (nSingles >= HowMany)
	printError("Too many ORF exons: decrease RSINGL parameter");
  
  return(nSingles);
}
