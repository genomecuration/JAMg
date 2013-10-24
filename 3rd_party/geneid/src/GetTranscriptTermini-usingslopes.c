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

/*  $Id: GetTranscriptTermini-usingslopes.c,v 1.3 2011/01/13 11:06:16 talioto Exp $  */

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
  long ns;
  float thresh = 1.2;
  /* Final number of potential transcript termini */
  ns = 0;  
/*   char mess[MAXSTRING]; */
  float score;
  long is;
  long j=0;
  long js=j;
  int window = 3;
  int ssExists = 0;
  int cluster_edge = 0;
  float pes=0;
  float splicethresh = 5;
  float prev_score =0;
  float prev_prev_score = 0;
  is = 0;     

  while ((is < (l2- l1 + 1)) && (ns<NUMSITES))
    { 
      /* is = 0..right */
      ssExists =0;
      score=0.0;
      pes = (PeakEdgeScore(l1 + is,Strand,external,l1,l2,6));
      cluster_edge = ClusterEdge(l1 + is,Strand,external,l1,l2);
      score = score + pes; 
/*       sprintf(mess,"pos %ld:%f",l1 + is -1,prev_score); */
/* 	  	printMess(mess); */
      if ((ns<NUMSITES)&&((cluster_edge == 1)||(prev_score>thresh && prev_score>prev_prev_score && prev_score>score))){
	j=js;
	while ((j < nAcceptors) && ((Acceptors+j)->Position < (l1 + is - 1 - window)))
	  j++;
	js=j;
	
	while ((j < nAcceptors) && ((Acceptors+j)->Position < (l1 + is - 1 + window))){
	  
		/* Save counter j for the next iteration */
	  if ((Acceptors+j)->Score > splicethresh){
	    ssExists++;
	  }
	  j++;
	  
	}
	if (ssExists == 0){
/* 	  sprintf(mess,"PES TSS: pos %ld:%f",LengthSequence - (l1 + is),score); */
/* 	  	printMess(mess); */
/* 	  sprintf(mess,"PES TSS: pos %ld:%f",l1 + is -1,prev_score); */
/* 	  	printMess(mess); */
	  sc[ns].Position=(l1 + is -1);  
	  sc[ns].Score=prev_score;
	  sc[ns].class=U2;
	  ns++;
	}else{
	  
/* 	  sprintf(mess,"SS exists PES TSS: pos %ld:%f",l1 + is -1,prev_score); */
/* 	  	printMess(mess); */
	}
      }
      is++;
      prev_prev_score = prev_score;
      prev_score = score;
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
  /* Final number of potential transcript termini */
  /* ns = 0;   */
/*   char mess[MAXSTRING]; */
  float score;
  long is;
  long j=0;
  long js=j;
  int window = 3;
  int ssExists = 0;
  int cluster_edge = 0;
  float thresh = 1.2;
  float pes=0;
  float splicethresh = 5;
  is = 0;     
  float prev_score =0;
  float prev_prev_score = 0;
  while ((is < (l2- l1 + 1)) && (ns<NUMSITES))
    { 
      /* is = 0..right */
      ssExists = 0;
      score=0.0;
      pes = (PeakEdgeScore(l1 + is,Strand,external,l1,l2,6));
      cluster_edge = ClusterEdge(l1 + is,Strand,external,l1,l2);
      score = score - pes; /* (pes>0?pes:-pes); */
		  
      if ((ns<NUMSITES)&&((cluster_edge == -1)||(prev_score>thresh && prev_score>prev_prev_score && prev_score>score))){
	j=js;
	while ((j < nDonors) && ((Donors+j)->Position < (l1 + is -1 - window)))
	  j++;
	js=j;
	
	while ((j < nDonors) && ((Donors+j)->Position < (l1 + is -1 + window))){
		/* Save counter j for the next iteration */
	  if ((Donors+j)->Score > splicethresh){
	    ssExists++;
	  }
	  j++;
	}
	if (ssExists == 0){
/* 	  sprintf(mess,"TES pos: %ld   PES TES: pos %ld:%f",l1 + is,LengthSequence - (l1 + is),score); */
/* 	  	printMess(mess); */
	  sc[ns].Position=(l1 + is -1);  
	  sc[ns].Score=prev_score;
	  sc[ns].class=U2;
	  ns++;
	}
      }
      prev_prev_score = prev_score;
      prev_score = score;
      is++;
    }

  if (ns >= NUMSITES)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);


}




