/*************************************************************************
*                                                                        *
*   Module: GetTranscriptTermini                                         *
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

/*  $Id: GetTranscriptTermini.c,v 1.2 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic sites */
extern long NUMSITES;

long GetTSS(
	    site* sc,
	    site* Acceptors, long nAcceptors,
	    packExternalInformation* external,
	    packHSP* hsp,
	    int Strand,
	    long LengthSequence, 
	    long l1,
	    long l2
	    ) 
{
  short x;
  short frameStart;
  long i = 0;
  long ns;
  long rl1;
  long rl2;
  long j=0;
  long js=j;
  int window = 10;
  long maxP2 = 0;
  /* Final number of potential transcript termini */
  ns = 0;  
/*   char mess[MAXSTRING]; */
  long tempcoord;

  x=0;
  if (Strand == FORWARD)
    {

    
/*   if (Strand == FORWARD) */
/*     { */
 
	  
      if (hsp != NULL)
	{
/* 	  i = external->iSegments[x]; */
/* 	  for (i = 0;  */
/* 	       i < hsp->nSegments[x]; i++){ */
/* /\* 	    sprintf(mess,"LengthSequence: %ld l1: %ld l2: %ld i: %ld p: %ld",LengthSequence,l1,l2,i, hsp->sPairs[x][i]->Pos1); *\/ */
/* /\* 	    printMess(mess); *\/ */
/* 	  } */
	  /* A. Skip HSPs out of this range: [l1,l2] */
	  for (i = 0;
/* 	       i = external->iSegments[x];  */
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1; 
	       i++)
	    ;
			  
	  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
	  for (
/* 	       i = external->iSegments[x] */
		 ; 
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1; 
	       i++)
	    {
	      /* Do nothing */
	    }
			  
	  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
	  for (; 
	       /* i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP; */ 
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2; 
	       i++)
	    {
	      if (hsp->sPairs[x][i]->Pos1>=l1 &&((i == 0) /* external->iSegments[x]) */
		  ||
						 ((maxP2 + UTRMAXGAP) < hsp->sPairs[x][i]->Pos1))){
		j=js;
		while ((j < nAcceptors) && ((Acceptors+j)->Position < (hsp->sPairs[x][i]->Pos1 - window)))
		  j++;
		
		/* Save counter j for the next iteration */
		js=j;
		if ((Acceptors+j)->Score < -1 || (Acceptors+j)->Position>=(hsp->sPairs[x][i]->Pos1 + window)){
		  if (ns<NUMSITES){
		    /* 		  sprintf(mess,"hsp->sPairs[x][i]->Pos1: %ld",hsp->sPairs[x][i]->Pos1); */
		    /* 		  printMess(mess); */
		    sc[ns].Position=hsp->sPairs[x][i]->Pos1 - COFFSET;  
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{
		  /* if ((Acceptors+j)->Position<(hsp->sPairs[x][i]->Pos1 + window)){ */
/* 		    (Acceptors+j)->Score = (Acceptors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos2 variable*/
	      maxP2 = MAX(hsp->sPairs[x][i]->Pos2,maxP2);
	    }
			  
	  /* Update partial counter: previous HSPs are useless for next split */
	  /* external->iSegments[x] = i; */ 	  
			  
	  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
	  for (; 
	       i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2; 
	       i++)
	    {
	      if ((i == external->iSegments[x])||((maxP2 + UTRMAXGAP) < hsp->sPairs[x][i]->Pos1)){
		j=js;
		while ((j < nAcceptors) && ((Acceptors+j)->Position < (hsp->sPairs[x][i]->Pos1 - window)))
		  j++;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Acceptors+j)->Score < -1 || (Acceptors+j)->Position>=(hsp->sPairs[x][i]->Pos1 + window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos1 - COFFSET;  
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{

/* 		  if ((Acceptors+j)->Position<(hsp->sPairs[x][i]->Pos1 + window)){ */
/* 		    (Acceptors+j)->Score = (Acceptors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos2 variable*/
	      maxP2 = MAX(hsp->sPairs[x][i]->Pos2,maxP2);
	    }
	}
		
    }
  else
    {
      /* HSPs in REVERSE strand */
      frameStart = FRAMES;
      tempcoord = l2;
      rl2 = LengthSequence - l1 + 1;
      rl1 = LengthSequence - l2 + 1;
      x=frameStart;
      js = nAcceptors;
      maxP2 = 0;
      /* HSPs in reverse strand are reverse-sorted by Position2 */
      
      if (hsp != NULL)
	{
/* 	  i = external->iSegments[x]; */
	  i = 0;
/* 	  for (i = 0; */
/* /\* 	       i = external->iSegments[x];  *\/ */
/* 	       i < hsp->nSegments[x]; i++){ */
/* /\* 	    sprintf(mess,"LengthSequence: %ld l1: %ld l2: %ld i: %ld p: %ld",LengthSequence,l1,l2,i, hsp->sPairs[x][i]->Pos1); *\/ */
/* /\* 	    printMess(mess); *\/ */
/* 	  } */
	  /* A. Skip HSPs out of this range: [l1,l2] */
	  for (i = hsp->nSegments[x] -1;
/* 	       i = external->iSegments[x]; */
	       i >= 0 && hsp->sPairs[x][i]->Pos2 < l1;
	       i--)
	    ;
			  
	  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
	  for (
/* 	       i = external->iSegments[x] */
		 ;
	       i >= 0 && hsp->sPairs[x][i]->Pos1 < l1;
	       i--)
	    {
	      if (hsp->sPairs[x][i]->Pos2 <= l2 && ((i == 0)||((maxP2 + UTRMAXGAP) < hsp->sPairs[x][i]->Pos1))){
		j=js;
		while ((j >= 0) && ((Acceptors+j)->Position > (hsp->sPairs[x][i]->Pos1 + window)))
		  j--;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Acceptors+j)->Score < -1 || (Acceptors+j)->Position<=(hsp->sPairs[x][i]->Pos1 - window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos1 - COFFSET;
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{
/* 		  sprintf(mess,"hsp->sPairs[x][i]->Pos1: %ld",hsp->sPairs[x][i]->Pos1); */
/* 		  printMess(mess); */
/* 		  if ((Acceptors+j)->Position>(hsp->sPairs[x][i]->Pos1 - window)){ */
/* 		    (Acceptors+j)->Score = (Acceptors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos2 variable*/
	      maxP2 = MAX(hsp->sPairs[x][i]->Pos2,maxP2);
	    }
			  
	  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
	  for (;
	       /* i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP; */
	       i >=0 && hsp->sPairs[x][i]->Pos2 <= l2;
	       i--)
	    {
/* 	      sprintf(mess,"hsp->sPairs[x][i]->Pos1: %ld",hsp->sPairs[x][i]->Pos1); */
/* 	      printMess(mess); */
	      if ((i == 0)||((maxP2 + UTRMAXGAP) < hsp->sPairs[x][i]->Pos1)){
		j=js;
		while ((j >= 0) && ((Acceptors+j)->Position > (hsp->sPairs[x][i]->Pos1 + window)))
		  j--;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Acceptors+j)->Score < -1 || (Acceptors+j)->Position<=(hsp->sPairs[x][i]->Pos1 - window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos1 - COFFSET;
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{
/* 		  sprintf(mess,"hsp->sPairs[x][i]->Pos1: %ld",hsp->sPairs[x][i]->Pos1); */
/* 		  printMess(mess); */
/* 		  if ((Acceptors+j)->Position>(hsp->sPairs[x][i]->Pos1 - window)){ */
/* 		    (Acceptors+j)->Score = (Acceptors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos2 variable*/
	      maxP2 = MAX(hsp->sPairs[x][i]->Pos2,maxP2);
	    }
			  
	  /* Update partial counter: previous HSPs are useless for next split */
/* 	  external->iSegments[x] = i; */
			  
	  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
	  for (;
	       i >=0 &&  hsp->sPairs[x][i]->Pos1 <= l2;
	       i--)
	    {
				  
	    }
	}
		
    }

  if (ns >= NUMSITES)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);


}

long GetTES(
	    site* sc,
	    site* Donors, long nDonors,
	    packExternalInformation* external,
	    packHSP* hsp,
	    int Strand,
	    long LengthSequence, 
	    long l1,
	    long l2,
	    long ns
	    ) 
{
  short x;
  short frameStart;
  long i = 0;
  /* long ns; */
  /* Final number of potential transcript termini */
  /* ns = 0; */  
  long tempcoord;
  long rl1;
  long rl2;
/*   char mess[MAXSTRING]; */
  long j = 0;
  long js = 0;
  int window = 10;
  long minP1 = l2;
  long maxP1 = 0;
  if (Strand == FORWARD)
    {

      x=0;
	  
      if (hsp != NULL)
	{
/* 	  i = external->iSegments[x]; */
/* 	  for (i = 0;  */
/* 	       i < hsp->nSegments[x]; i++){ */
/* /\* 	    sprintf(mess,"LengthSequence: %ld l1: %ld l2: %ld i: %ld p: %ld",LengthSequence,l1,l2,i, hsp->sPairs[x][i]->Pos1); *\/ */
/* /\* 	    printMess(mess); *\/ */
/* 	  } */
	  /* A. Skip HSPs out of this range: [l1,l2] */
	  for (i = hsp->nSegments[x]-1; 
	       i >= 0 && hsp->sPairs[x][i]->Pos1 > l2; 
	       i--)
	    ;
			  
	  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
	  for (; 
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2; 
	       i--)
	    {
	      if (hsp->sPairs[x][i]->Pos1 <= l2 && ((i == hsp->nSegments[x]-1)||((minP1 - UTRMAXGAP) > hsp->sPairs[x][i]->Pos2))){
		j=js;
		while ((j < nDonors) && ((Donors+j)->Position <= (hsp->sPairs[x][i]->Pos1 - window)))
		  j++;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Donors+j)->Score < -1 || (Donors+j)->Position>(hsp->sPairs[x][i]->Pos1 + window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos2 - COFFSET;  
		    sc[ns].Score=0;
		    sc[ns].class=U2;

/* 		    sprintf(mess,"TES pos: %ld  l1:%ld",sc[ns].Position,l1); */
/* 		    printMess(mess); */
		    ns++;
		  }
		}else{
/* 		  if ((Donors+j)->Position<=(hsp->sPairs[x][i]->Pos1 + window)){ */
/* 		    (Donors+j)->Score = (Donors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos1 variable*/
	      minP1 = MIN(hsp->sPairs[x][i]->Pos1,minP1);
	    }
			  
	  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
	  for (; 
	       i >= 0 && hsp->sPairs[x][i]->Pos1 >= l1; 
	       i--)
	    {
	      if ((i == hsp->nSegments[x]-1)||((minP1 - UTRMAXGAP) > hsp->sPairs[x][i]->Pos2)){
		j=js;
		while ((j < nDonors) && ((Donors+j)->Position <= (hsp->sPairs[x][i]->Pos1 - window)))
		  j++;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Donors+j)->Score < -1 || (Donors+j)->Position>(hsp->sPairs[x][i]->Pos1 + window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos2 - COFFSET;  
		    sc[ns].Score=0;
		    sc[ns].class=U2;
/* 		    sprintf(mess,"TES pos (pos1>l1): %ld  l1:%ld",sc[ns].Position,l1); */
/* 		    printMess(mess); */
		    ns++;
		  }
		}else{
/* 		  if ((Donors+j)->Position<=(hsp->sPairs[x][i]->Pos1 + window)){ */
/* 		    (Donors+j)->Score = (Donors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos1 variable*/
	      minP1 = MIN(hsp->sPairs[x][i]->Pos1,minP1);
	    }
			  
	  /* Update partial counter: previous HSPs are useless for next split */
	  /* external->iSegments[x] = i;  */	  
			  
	  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
	  for (; 
	       i >=0 &&  hsp->sPairs[x][i]->Pos2 >= l1; 
	       i--)
	    {

				  
	    }
	}
		
    }
  else
    {
      /* HSPs in REVERSE strand */
 
      frameStart = FRAMES;
      tempcoord = l2;
      rl2 = LengthSequence - l1 + 1;
      rl1 = LengthSequence - l2 + 1;
      x=frameStart;
      js = nDonors;
      maxP1 =0;
      /* HSPs in reverse strand are reverse-sorted by Position2 */

      if (hsp != NULL)
	{
/* 	  i = external->iSegments[x]; */
	  i = 0;
	 /*  for (i = 0;  */
/* 	       i < hsp->nSegments[x]; i++){ */
/* /\* 	    sprintf(mess,"LengthSequence: %ld l1: %ld l2: %ld i: %ld p: %ld",LengthSequence,l1,l2,i, hsp->sPairs[x][i]->Pos1); *\/ */
/* /\* 	    printMess(mess); *\/ */
/* 	  } */
	  /* A. Skip HSPs out of this range: [l1,l2] */
	  for (i = 0;
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2;
	       i++)
	    ;
			  
	  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
	  for (
/* 	       i = external->iSegments[x] */
		 ;
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2;
	       i++)
	    {

	    }
			  
	  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
	  for (;
	       i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1;
	       i++)
	    {
	      if (hsp->sPairs[x][i]->Pos2 <= l2 &&((i == 0)||((maxP1 - UTRMAXGAP) > hsp->sPairs[x][i]->Pos2))){
		j=js;
		while ((j >= 0) && ((Donors+j)->Position >= (hsp->sPairs[x][i]->Pos1 + window)))
		  j--;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Donors+j)->Score < -1 || (Donors+j)->Position<(hsp->sPairs[x][i]->Pos1 - window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos2 - COFFSET;
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{
/* 		  if ((Donors+j)->Position>=(hsp->sPairs[x][i]->Pos1 - window)){ */
/* 		    (Donors+j)->Score = (Donors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos1 variable*/
	      maxP1 = MAX(hsp->sPairs[x][i]->Pos1,maxP1);
	    }
			  
	  /* Update partial counter: previous HSPs are useless for next split */
/* 	  external->iSegments[x] = i; */
			  
	  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
	  for (;
	       i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1;
	       i++)
	    {
	      if (hsp->sPairs[x][i]->Pos2 <= l2 && ((i == external->iSegments[x])||((maxP1 - UTRMAXGAP) > hsp->sPairs[x][i]->Pos2))){
		j=js;
		while ((j >= 0) && ((Donors+j)->Position >= (hsp->sPairs[x][i]->Pos1 + window)))
		  j--;
	  
		/* Save counter j for the next iteration */
		js=j;
		if ((Donors+j)->Score < -1 || (Donors+j)->Position<(hsp->sPairs[x][i]->Pos1 - window)){
		  if (ns<NUMSITES){
		    sc[ns].Position=hsp->sPairs[x][i]->Pos2 - COFFSET;
		    sc[ns].Score=0;
		    sc[ns].class=U2;
		    ns++;
		  }
		}else{
/* 		  if ((Donors+j)->Position>=(hsp->sPairs[x][i]->Pos1 - window)){ */
/* 		    (Donors+j)->Score = (Donors+j)->Score + 1; */
/* 		  } */
		}
	      }
	      /*update max pos1 variable*/
	      maxP1 = MAX(hsp->sPairs[x][i]->Pos1,maxP1);
	    }
	}
		
    }

  if (ns >= NUMSITES)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);


}




