/*************************************************************************
 *                                                                        *
 *   Module: BuildUTRExons                                                *
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

/*  $Id: BuildUTRExons.c,v 1.3 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (divided by RFIRST) */
/* Sequence is used to save information to prevent Stop codons in frame */
/* extern long NUMEXONS; */
extern long MAXBACKUPSITES;

long BuildUTRExons(
				site *Start, long nStarts, 
				site *Donor, long nDonors,
				int MaxDonors,
				int MaxExonLength,long cutPoint,
				char* ExonType,
				exonGFF *Exon,long nexons) 
{
  /* Best exons built by using the current start codon */
  exonGFF *LocalExon;
  int nLocalExons, LowestLocalExon;
  float LowestLocalScore;
/*   char mess[MAXSTRING]; */

  /* Maximum allowed number of predicted initial exons per fragment */
  long HowMany;
  
/*   int Frame; */
/*   int js; */

  long i, j, k, ks;
  int l;
  float pen = 0.002;
  /* Final number of predicted UTR exons */
  long nExon;
  
  /* Allocating space for exons built by using the current start codon */
  /* MaxDonors is the maximum allowed number of exons with this signal */
  if ((LocalExon = (exonGFF*) calloc(MaxDonors,sizeof(exonGFF))) == NULL)
    printError("Not enough memory: local UTR exons"); 
  
  /* Main loop, forall beginning sites looking for ending sites... */
  HowMany = (MAXBACKUPSITES)? (long)(nexons/RUTR): nexons;
  for (i = 0, j = 0, k = 0, nExon = 0;
       (nExon < HowMany) && (i < nStarts); i++)
    { 
      /* Reset the best local exons array */
      nLocalExons = 0;
      LowestLocalScore = INF;
      LowestLocalExon = 0;
	  
      /* Skip previous Donors to current Start */
      while (((Donor+k)->Position < (Start+i)->Position + 1) && (k < nDonors)){
/* 	sprintf(mess,"start pos: %ld\nend pos: %ld",(Start+i)->Position,(Donor+k)->Position); */
/* 	printMess(mess); */
	k++;
      }
	  
      /* Save counter k for the next iteration */
      ks=k;
	  
      /* a) Every Donor between Start and that Stop defines an initial exon */
      /* b) If not any Stop in frame after current Start: every donor is OK */
      while (((Donor+ks)->Position < (Start+i)->Position + MaxExonLength)
	     && (ks < nDonors)  && ((Donor+ks)->Position > (Start+i)->Position)
/* 	     &&((!strcmp(ExonType,sUTRTERMINAL)||!strcmp(ExonType,sUTRTERMINALHALF))?(Donor+ks)->Position >= cutPoint:1) */
	  
	     ) /* && ((Donor+ks)->Position > (Start+i)->Position+2)*/
	{
	  /* a. There is room to save this new exon */
	  if (nLocalExons < MaxDonors) 
	    {
	      
/* 	      if (!strcmp(ExonType,sUTRTERMINALHALF)){ */
/* 		sprintf(mess,"start pos: %ld\nend pos: %ld",(Start+i)->Position,(Donor+ks)->Position);  	      	printMess(mess);  */
/* 		  } */
	      /* 	sprintf(mess,"start pos: %ld\nend pos: %ld",(Start+i)->Position,(Donor+ks)->Position); */
/* 	      	printMess(mess); */
	      (LocalExon+nLocalExons)->Acceptor=(Start+i);
	      (LocalExon+nLocalExons)->Donor=(Donor+ks);
			  
	      /* Updating the worst exon pointer */
	      if ((Donor+ks)->Score - pen*((Donor+ks)->Position - (Start+i)->Position + 1) < LowestLocalScore)
		{
		  LowestLocalScore = (Donor+ks)->Score - pen*((Donor+ks)->Position - (Start+i)->Position + 1);
		  LowestLocalExon = nLocalExons;
		}
/* 	      sprintf(mess,"ExonType: %s   start pos: %ld   don pos: %ld   curr exon don score: %f    lowest exon donor score: %f",ExonType, (Start+i)->Position,(Donor+ks)->Position,((Donor+ks)->Score - pen*((Donor+ks)->Position - (Start+i)->Position + 1)),LowestLocalScore); */
/* 	      printMess(mess); */
	      nLocalExons++;
	    }
	  else 
	    {
	      /* b. Keep only the top scoring MaxDonors */
	      
	      if (((Donor+ks)->Score - pen*((Donor+ks)->Position - (Start+i)->Position + 1)) > LowestLocalScore)
		{
/* 		  sprintf(mess,"ExonType: %s   start pos: %ld   don pos: %ld   curr exon don score: %f    lowest exon donor score: %f",ExonType, (Start+i)->Position,(Donor+ks)->Position,((Donor+ks)->Score - pen*((Donor+ks)->Position - (Start+i)->Position + 1)),LowestLocalScore); */
/* 	      printMess(mess); */
		  /* Extract the worst exon and input the new exon */
		  for (l=LowestLocalExon;l<nLocalExons-1;l++)
		    LocalExon[l]=LocalExon[l+1];
				  
		  (LocalExon+nLocalExons-1)->Acceptor=(Start+i);
		  (LocalExon+nLocalExons-1)->Donor=(Donor+ks);
				  
		  /* Updating the worst exon pointer */
		  LowestLocalExon = 0;
		  LowestLocalScore = (LocalExon+0)->Donor->Score - pen*((LocalExon+0)->Donor->Position - (LocalExon+0)->Acceptor->Position + 1);
		  for (l=1;l<nLocalExons;l++)
		    if ((LocalExon+l)->Donor->Score - pen*((LocalExon+1)->Donor->Position - (LocalExon+1)->Acceptor->Position + 1)< LowestLocalScore) 
		      {
			LowestLocalScore = (LocalExon+l)->Donor->Score - pen*((LocalExon+1)->Donor->Position - (LocalExon+1)->Acceptor->Position + 1);
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
	  (Exon+nExon)->Remainder = 0; 
	  /* Assign Type according to donor subtype */
	  strcpy((Exon+nExon)->Type,ExonType);
	  strcpy((Exon+nExon)->Group,NOGROUP);
	  (Exon+nExon)->evidence = 0;		    
	  (Exon+nExon)->rValue = 0;
	  (Exon+nExon)->lValue = 0;
/* 	  sprintf(mess,"ExonType: %s   start pos: %ld   don pos: %ld   curr exon don score: %f",ExonType, (LocalExon+l)->Acceptor->Position,(LocalExon+l)->Donor->Position,(LocalExon+l)->Donor->Score); */
/* 	  printMess(mess); */
	  nExon++;
	  
	}
    }
  
  if (nExon >= HowMany)
    printError("Too many UTR exons: decrease RUTR parameter");
  
  free(LocalExon);
  
  return(nExon);
}

