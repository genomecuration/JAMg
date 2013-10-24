/*************************************************************************
*                                                                        *
*   Module: BuildInitialExons                                            *
*                                                                        *
*   From start/stop codons and donor sites, to build initial exons       *
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

/*  $Id: BuildInitialExons.c,v 1.10 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (divided by RFIRST) */
/* Sequence is used to save information to prevent Stop codons in frame */
/* extern long NUMEXONS; */
extern long MAXBACKUPSITES;

long BuildInitialExons(site *Start, long nStarts, 
                       site *Donor, long nDonors,
                       site *Stop, long nStops,
                       int MaxDonors,
					   char* ExonType,
					   char* Sequence,
                       exonGFF *Exon,long nexons) 
{
  /* Best exons built by using the current start codon */
  exonGFF *LocalExon;
  int nLocalExons, LowestLocalExon;
  float LowestLocalScore;
  
  /* Maximum allowed number of predicted initial exons per fragment */
  long HowMany;
  
  int Frame;
  long i, j, js, k, ks;
  int l;
  
  /* Final number of predicted initial exons */
  long nExon;
  
  
  /* Allocating space for exons built by using the current start codon */
  /* MaxDonors is the maximum allowed number of exons with this signal */
  if ((LocalExon = (exonGFF*) calloc(MaxDonors,sizeof(exonGFF))) == NULL)
    printError("Not enough memory: local first exons"); 
  
  /* Main loop, forall start codon looking for donor sites... */
  /* ...until the first stop codon in frame is reached */
  HowMany = (MAXBACKUPSITES)? (long)(nexons/RFIRST): nexons;
  for (i = 0, j = 0, k = 0, nExon = 0;
       (nExon < HowMany) && (i < nStarts); i++)
	{ 
	  /* Reset the best local exons array */
	  nLocalExons = 0;
	  LowestLocalScore = INF;
	  LowestLocalExon = 0;
	  
	  /* Computing the frame to look for the first stop in frame */
	  Frame = (Start+i)->Position % 3;
	  
	  /* Skip previous Stops to current Start */
	  while (((Stop+j)->Position+1 < (Start+i)->Position) && (j < nStops))
		j++;
	  
	  /* Save counter j for the next iteration */
	  js=j;
	  
	  /* Finding first Stop in Frame with current Start */
	  while ((((Stop+js)->Position+1) % 3 != Frame) && (js < nStops))
		js++;
	  
	  /* Skip previous Donors to current Start */
	  while (((Donor+k)->Position < (Start+i)->Position+2) && (k < nDonors))
		k++;
	  
	  /* Save counter k for the next iteration */
	  ks=k;
	  
	  /* a) Every Donor between Start and that Stop defines an initial exon */
	  /* b) If not any Stop in frame after current Start: every donor is OK */
	  while ((js == nStops || (Donor+ks)->Position < (Stop+js)->Position + 1 + 2)
		 && (ks < nDonors)) /* && ((Donor+ks)->Position > (Start+i)->Position+2)*/
		{
		  /* a. There is room to save this new exon */
		  if (nLocalExons < MaxDonors) 
			{
			  (LocalExon+nLocalExons)->Acceptor=(Start+i);
			  (LocalExon+nLocalExons)->Donor=(Donor+ks);
			  
			  /* Updating the worst exon pointer */
			  if ((Donor+ks)->Score < LowestLocalScore)
				{
				  LowestLocalScore = (Donor+ks)->Score;
				  LowestLocalExon = nLocalExons;
				}
			  nLocalExons++;
			}
		  else 
			{
			  /* b. Keep only the top scoring MaxDonors */
			  if ((Donor+ks)->Score > LowestLocalScore)
				{
				  /* Extract the worst exon and input the new exon */
				  for (l=LowestLocalExon;l<nLocalExons-1;l++)
					LocalExon[l]=LocalExon[l+1];
				  
				  (LocalExon+nLocalExons-1)->Acceptor=(Start+i);
				  (LocalExon+nLocalExons-1)->Donor=(Donor+ks);
				  
				  /* Updating the worst exon pointer */
				  LowestLocalExon = 0;
				  LowestLocalScore = (LocalExon+0)->Donor->Score;
				  for (l=1;l<nLocalExons;l++)
					if ((LocalExon+l)->Donor->Score < LowestLocalScore) 
					  {
						LowestLocalScore = (LocalExon+l)->Donor->Score;
						LowestLocalExon = l;
					  }
				} /* end if */
			} /* end else */
		  
		  /* Next donor between current start and stop in frame */
		  ks++;
		} /* end while */
	  
	  /* Save the best exons beginning by the current start */
	  for (l=0;(l<nLocalExons) && (nExon<HowMany);l++) 
		{
		  Exon[nExon] = LocalExon[l];
		  (Exon+nExon)->Frame = 0;
		  (Exon+nExon)->Remainder = ((Exon+nExon)->Donor->Position -
									 (Exon+nExon)->Acceptor->Position + 1 ) % 3;
		  (Exon+nExon)->Remainder = (3 - (Exon+nExon)->Remainder) % 3; 
         /* Assign Type according to donor subtype */
		 strcpy((Exon+nExon)->Type,ExonType);
         strcpy((Exon+nExon)->Group,NOGROUP);
         (Exon+nExon)->evidence = 0;

		 /* Store info to prevent building stops in frame */
		 ComputeStopInfo((Exon+nExon),Sequence);

         nExon++;
		}
	}
  
  if (nExon >= HowMany)
	printError("Too many initial exons: decrease RFIRST parameter");
  
  free(LocalExon);
  
  return(nExon);
}

