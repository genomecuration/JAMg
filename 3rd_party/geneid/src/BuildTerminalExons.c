/*************************************************************************
*                                                                        *
*   Module: BuildTerminalExons                                           *
*                                                                        *
*   From acceptor sites and stop codons, to build initial exons          *
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

/*  $Id: BuildTerminalExons.c,v 1.9 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (divided by RTERMI) */
/* Sequence is used to save information to prevent Stop codons in frame */
/* extern long NUMEXONS; */
extern long MAXBACKUPSITES;

long BuildTerminalExons (site *Acceptor, long nAcceptors, 
                         site *Stop, long nStops,
                         long LengthSequence,
                         long cutPoint,
						 char* ExonType,
						 char* Sequence,
						 exonGFF* Exon, long nexons)
{
  int Frame[FRAMES];
  long i, f, j, js;
/*   char mess[MAXSTRING]; */
  /* Final number of predicted terminal exons */
  long nExon;
  
  /* Maximum allowed number of predicted initial exons per fragment */
  long HowMany;
  
  
  /* Main loop: forall Acceptor, first Stop in every Frame defines an exon */
  /* There are therefore, 3 terminal exons starting by this Acceptor */
  HowMany = (MAXBACKUPSITES)? (long)(nexons/RTERMI): nexons;
  
  for (i=0, j=0, js=0, nExon=0;
       (i<nAcceptors) && (nExon<HowMany);
       i++)
    {
      /* Open the 3 frames for current Acceptor */
      for (f=0;f<FRAMES;f++)
		Frame[f]=1;
	  
      /* Skip previous Stops to current Acceptor */
      while (((Stop+j)->Position+1 < (Acceptor+i)->Position) && (j < nStops))
		j++;
	  
      /* Save counter j for the next iteration */
      js=j;
      
      /* Use current Stops if its frame is still opened */
      while ((Frame[0]==1 || Frame[1]==1 || Frame[2]==1)
			 && (js < nStops)
			 && (nExon<HowMany))
		{
		  if (Frame[f=((Stop+js)->Position - (Acceptor+i)->Position + 1) % 3])
			{
			  /* Save the new exon for the current couple (Acceptor,Stop) */
			  /* CutPoint: to preserve sorted exons between fragments */
			  if ((Stop+js)->Position >= (Acceptor+i)->Position 
			      &&(Stop+js)->Position >= cutPoint
			      )
				{
				  
				  (Exon+nExon)->Acceptor=(Acceptor+i);
				  (Exon+nExon)->Donor=(Stop+js);
				  (Exon+nExon)->Frame = f;
				  (Exon+nExon)->Remainder = 0; 
				  /* Assign Type according to acceptor subtype */
				  strcpy((Exon+nExon)->Type,ExonType);
				  strcpy((Exon+nExon)->Group,NOGROUP);
				  (Exon+nExon)->evidence = 0;

				  /* Store info about frame and remainder nucleotides to avoid building stops in frame */
				  ComputeStopInfo((Exon+nExon),Sequence);
				  nExon++;
				}
			  /* Close this frame */
			  Frame[f]=0;
			} 
		  js++;
		} /* next stop */     
    } /* next acceptor */   
  
  if (nExon >= HowMany)
	printError("Too many terminal exons: decrease RTERMI parameter");
  
  return(nExon);
}
