/*************************************************************************
*                                                                        *
*   Module: BuildInternalExons                                           *
*                                                                        *
*   From acceptor/donor sites and stop codons, to build internal exons   *
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

/*  $Id: BuildZeroLengthExons.c,v 1.2 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Sequence is used to save information to prevent Stop codons in frame */
/* Maximum allowed number of generic exons (divided by RINTER) */
/* extern long NUMEXONS; */
extern long MAXBACKUPSITES;
extern int RSS;
extern float RSSDON;
extern float RSSACC;

long BuildZeroLengthExons(site *Acceptor, long nAcceptors, 
                        site *Donor, long nDonors,
                        site *Stop, long nStops,
                        int MaxDonors,
			char* ExonType,
			char* Sequence,
                        exonGFF* Exon, long nexons)
{
  /* Best exons built using the current Acceptor and frame availability */
  struct iexon
  {
    site *Acceptor;
    site *Donor;
    int Frame[FRAMES];
  } *LocalExon; 
  int nLocalExons, LowestLocalExon;
  float LowestLocalScore;
  /* char mess[MAXSTRING]; */

  /* Boolean array of windows: closed or opened */ 
  int Frame[FRAMES];
  
  long i, j, js, k;
  int f, l, ll;
  
  /* Maximum allowed number of predicted internal exons per fragment */
  long HowMany;
  
  /* Final number of predicted internal exons */
  long nExon;
  
  /* Allocating space for exons built by using the current Acceptor */
  /* MaxDonors is the maximum allowed number of exons with this signal */
  if ((LocalExon =
       (struct iexon *) calloc(MaxDonors,sizeof(struct iexon))) == NULL)
    printError("Not enough memory: local initial exons"); 
  
  /* Main loop, forall Acceptor looking for donor sites... */
  /* ...until 3 windows are closed due to 3 stop codons */
  HowMany = (MAXBACKUPSITES)? (long)(nexons/RINTER): nexons;
  
  
  for (i=0, j=0, k=0, nExon = 0;
       (i < nAcceptors) && (nExon<HowMany); i++)
    {
      /* Open the 3 windows */
      for (f=0;f<FRAMES;f++) 
	Frame[f]=1;

      /* Reset the best local exons array */
      nLocalExons = 0;
      LowestLocalScore = INF;
      LowestLocalExon = 0;

      /* Skip previous Stops to current Acceptor */
      while ((j < nStops) && ((Stop+j)->Position+1 < (Acceptor+i)->Position))
	j++;
	  
      /* Save counter j for the next iteration */
      js=j;
      
      /* Skip previous Donors to current Acceptor */
      while ((k < nDonors) && ((Donor+k)->Position < (Acceptor+i)->Position - 1))
	k++;
      if (RSS){
	if (((Donor+k)->Position == (Acceptor+i)->Position - 1)&&((Donor+k)->Score > RSSDON) &&((Acceptor+i)->Score > RSSACC)){
	  /* Make a zero length exon representing the recursive splice site */
	  (LocalExon+nLocalExons)->Acceptor=(Acceptor+i);
	  (LocalExon+nLocalExons)->Donor=(Donor+k);
					  
	  /* Saving the exon in the opened frames */
	  for (ll=0;ll<FRAMES;ll++)
	    (LocalExon+nLocalExons)->Frame[ll]=Frame[ll];
					  
	  /* Updating the worst exon pointer */
	  if ((Donor+k)->Score < LowestLocalScore)
	    {
	      LowestLocalScore = (Donor+k)->Score;
	      LowestLocalExon = nLocalExons;
	    }
	  nLocalExons++;
	}
      }
     	  	  
      /* Save predicted exons for the current Acceptor site */
      for (l=0;(l<nLocalExons) && (nExon<HowMany);l++) 
	{
	  /* There are exons in every frame */
	  for (ll=0;(ll<FRAMES) && (nExon < HowMany);ll++)
	    if ((LocalExon+l)->Frame[ll]) 
	      {
		/* Frame was opened then */
		(Exon+nExon)->Acceptor = (LocalExon+l)->Acceptor;
		(Exon+nExon)->Donor = (LocalExon+l)->Donor;
		(Exon+nExon)->Frame = ll;
		(Exon+nExon)->Remainder = ((Exon+nExon)->Donor->Position -
					   ((Exon+nExon)->Acceptor->Position +
					    (Exon+nExon)->Frame)+ 1 ) % 3;
		(Exon+nExon)->Remainder = (3 - (Exon+nExon)->Remainder) % 3;
		strcpy((Exon+nExon)->Type,ExonType);
		strcpy((Exon+nExon)->Group,NOGROUP);
		(Exon+nExon)->evidence = 0;
				
		/* Save info (frame/remainder bases to avoid stops in frame */
		if(!RSS){
		  ComputeStopInfo((Exon+nExon),Sequence);
		} else {
		  if ((Exon+nExon)->Donor->Position == (Exon+nExon)->Acceptor->Position - 1){
		    (Exon+nExon)->rValue = 0;
		    (Exon+nExon)->lValue = 0;
		  }else{
		    ComputeStopInfo((Exon+nExon),Sequence);
		  }
		}
		nExon++;
	      }
	}
    }
  
  if (nExon >= HowMany)
    printError("Too many zero-length exons: decrease RINTER parameter");
  
  free(LocalExon);
  
  return(nExon);
}


