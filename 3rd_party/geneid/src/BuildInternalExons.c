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

/*  $Id: BuildInternalExons.c,v 1.11 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Sequence is used to save information to prevent Stop codons in frame */
/* Maximum allowed number of generic exons (divided by RINTER) */
/* extern long NUMEXONS; */
extern long MAXBACKUPSITES;
extern int RSS;

long BuildInternalExons(site *Acceptor, long nAcceptors, 
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
  
  /* Boolean array of windows: closed or opened */ 
  int Frame[FRAMES];
  
  long i, j, js, k, ks;
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
      
      /* Skip previous Donors to current Acceptor plus minimum exon length */
      /* Minimal length: EXONLENGTH */
      while ((k < nDonors) && ((Donor+k)->Position < (Acceptor+i)->Position + EXONLENGTH))
	k++;
	  
      /* Save counter k for the next iteration */
      ks=k;
      
      /* A. While there are opened frames, look for next Stop to close them */
      /* Donors between Acceptor and current Stop define internal exons */
      while ((js < nStops) && (Frame[0]==1 || Frame[1]==1 || Frame[2]==1))
	{
	  /* If this window is still opened... */
	  if ((Frame[f=((Stop+js)->Position - (Acceptor+i)->Position + 1) % 3]))
	    {
	      /* Donors between Acceptor and Stop defines internal exons */
	      while (((Donor+ks)->Position < (Stop+js)->Position + 1 + 2)
		     && (ks < nDonors)) /*&&((Donor+ks)->Position > (Acceptor+i)->Position + EXONLENGTH)*/
		{
		  /* a. There is room to save this new exon */
		  if (nLocalExons < MaxDonors)
		    {
		      (LocalExon+nLocalExons)->Acceptor=(Acceptor+i);
		      (LocalExon+nLocalExons)->Donor=(Donor+ks);
					  
		      /* Saving the exon in the opened frames */
		      for (ll=0;ll<FRAMES;ll++)
			(LocalExon+nLocalExons)->Frame[ll]=Frame[ll];
					  
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
		      if ((Donor+ks)->Score >= LowestLocalScore)
			{
			  /* Extract the worst exon and input the new exon */
			  for (l=LowestLocalExon;l<nLocalExons-1;l++)
			    LocalExon[l]=LocalExon[l+1];
			  (LocalExon+nLocalExons-1)->Acceptor=(Acceptor+i);
			  (LocalExon+nLocalExons-1)->Donor=(Donor+ks);
						  
			  /* Saving the exon in the opened frames */
			  for (ll=0;ll<FRAMES;ll++)
			    (LocalExon+nLocalExons-1)->Frame[ll]=Frame[ll];
						  
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
		  /* Next donor */
		  ks++;
		} /* while there are donors left */
			  
	      /* Close the window corresponding to the Stop found */
	      Frame[f]=0;
	    } /* If this window was still opened... */
	  /* Next Stop */
	  js++;
	} /* end while-frames opened scan stops */
	  
	  /* B. EXTRA: Stop list finished but there are frames still opened... */
      if ((js == nStops) && (Frame[0]==1 || Frame[1]==1 || Frame[2]==1))
	{
	  /* if any frame is opened and there are not more stops left */
	  /* then every donor forms an internal exon with the Acceptor */
	  while ((ks < nDonors)) /*&&((Donor+ks)->Position > (Acceptor+i)->Position + EXONLENGTH)*/
	    {
	      /* a. There is room to save this new exon */
	      if (nLocalExons < MaxDonors)
		{
		  (LocalExon+nLocalExons)->Acceptor=(Acceptor+i);
		  (LocalExon+nLocalExons)->Donor=(Donor+ks);
				  
		  /* Saving the exon in the opened frames */
		  for (ll=0;ll<FRAMES;ll++)
		    (LocalExon+nLocalExons)->Frame[ll]=Frame[ll];
				  
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
		  if ((Donor+ks)->Score >= LowestLocalScore)
		    {
		      /* Extract the worst exon and input the new exon */
		      for (l=LowestLocalExon;l<nLocalExons-1;l++)
			LocalExon[l]=LocalExon[l+1];
		      (LocalExon+nLocalExons-1)->Acceptor=(Acceptor+i);
		      (LocalExon+nLocalExons-1)->Donor=(Donor+ks);
					  
		      /* Saving the exon in the opened frames */
		      for (ll=0;ll<FRAMES;ll++)
			(LocalExon+nLocalExons-1)->Frame[ll]=Frame[ll];
					  
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
		}/* end else */  
	      /* Next donor */
	      ks++;
	    } /* while there are donors... */
	}/* end if EXTRA */
	  
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
		  }
		}
		nExon++;
	      }
	}
    }
  
  if (nExon >= HowMany)
    printError("Too many internal exons: decrease RINTER parameter");
  
  free(LocalExon);
  
  return(nExon);
}


