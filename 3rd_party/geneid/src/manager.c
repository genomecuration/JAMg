/*************************************************************************
*                                                                        *
*   Module: manager                                                      *
*                                                                        *
*   Management of prediction actions: signals, exons and scores          *
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

/* $Id: manager.c,v 1.15 2011/01/13 11:06:16 talioto Exp $ */

#include "geneid.h"

extern int scanORF;
extern int PAS;
extern int U12GTAG;
extern int U12ATAC;
extern int U2GCAG;
extern int U2GTA;
extern int U2GTG;
extern int U2GTY;
extern int RSS;
extern int GENAMIC;
extern int EFP;
extern int EIP;
extern int ETP;
extern int ESP;
extern int EOP;
extern int EXP;
extern int SRP;
extern int UTR;

extern long NUMSITES,NUMEXONS;

/* Management of splice sites prediction and exon construction/scoring */
void  manager(char *Sequence, 
	      long LengthSequence,
	      packSites* allSites,
	      packExons* allExons,
	      long l1, long l2, long lowerlimit, long upperlimit,
	      int Strand,
	      packExternalInformation* external,
	      packHSP* hsp,
	      gparam* gp,
	      gparam** isochores,
	      int nIsochores,
	      packGC* GCInfo,
	      site* acceptorsites,
	      site* donorsites,
	      site* tssites,
	      site* tesites
	      )
{
  char mess[MAXSTRING];

  /* For sorting sites */
/*   site* acceptorsites;  */
/*   site* donorsites;  */
  long l1a, l1b,
	l2a, l2b,
	l1c, l2c;

  long cutPoint;
 
  /* 0. Define boundaries of splice site prediction
	 according to current split positions and strand selected */
  if (Strand == FORWARD)
    {
      /* Forward sense */
      /* Start codons and Acceptor sites limits */
      l1a = l1;
      l2a = (l2 == upperlimit)? l2 : l2 - OVERLAP;

      /* Donor sites limits */
      l1b = l1;
      l2b = l2;

      /* Stop codon limits */
      l1c = l1;
      l2c = l2;

      /* Terminal/Single exons: */
      /* are allowed if their Stop codon is placed behind cutPoint */
      /* FWD: every stop codon might be used without problems */
      cutPoint = l1;
    }
  else
    {
      /* Reverse sense */
      /* Start codons and Acceptor sites limits */
      l1a = l1;
      l2a = l2;

      /* Donor sites limits */
      l1b = (l1 == lowerlimit)? l1: l1 + OVERLAP;
      l2b = l2;

      /* Stop codon limits */
      l1c = l1;
      l2c = l2;

      /* Terminal/Single exons: */
      /* are allowed if their Stop codon is placed behind cutPoint (RVS) */
      /* RVS: reading from right to left the forward sense sequence */
      cutPoint = (l1 == lowerlimit)? l1 : l1 + OVERLAP;
    }

/*   sprintf(mess,"Strand:%i\nl1a:%ld\nl1b:%ld\nl2a:%ld\nl2b:%ld\nl1c:%ld\nl2c:%ld\ncutPoint:%ld\n",Strand,l1a, l1b,l2a, l2b,l1c, l2c,cutPoint); */
/*   printMess(mess); */

  /* 0. Preprocss HSPs */
  if (SRP){
    ProcessHSPs(l1, l2, Strand, 
		external, hsp);
  }

  /* 1. Predicting splice sites of current split of DNA sequence */ 
  printMess ("Computing sites ...");

  allSites->nStartCodons =
    GetSitesWithProfile(Sequence,gp->StartProfile,allSites->StartCodons,l1a,l2a);
  sprintf(mess, "Start Codons \t\t%8ld", allSites->nStartCodons);
  printRes(mess);
  
  long numAccsites = 0;
  
  allSites->nAcceptorSites =
    BuildAcceptors(Sequence,
		   U2,
		   sU2type,
		   sU2,
		   gp->AcceptorProfile,
		   gp->PolyPTractProfile,
		   gp->BranchPointProfile,
		   allSites->AcceptorSites,
		   l1a,l2a,numAccsites,NUMSITES,Strand,external);
  
  sprintf(mess, "Acceptor Sites \t\t%8ld", allSites->nAcceptorSites - numAccsites);
  numAccsites = allSites->nAcceptorSites;
  printRes(mess);

  if (U12GTAG){ 
	  allSites->nAcceptorSites =
	    BuildU12Acceptors(Sequence,U12gtag,sU12type,
					   sU12gtag,
					   gp->U12gtagAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->AcceptorSites,
					   l1a,l2a,numAccsites,NUMSITES,Strand,external);

	  sprintf(mess, "U12gtag Acceptor Sites \t%8ld", allSites->nAcceptorSites - numAccsites);
	  numAccsites = allSites->nAcceptorSites;
	  printRes(mess);
  }
  if (U12ATAC){ 
	  allSites->nAcceptorSites =
	    BuildU12Acceptors(Sequence,U12atac,sU12type,
				 	   sU12atac,
					   gp->U12atacAcceptorProfile,
					   gp->U12BranchPointProfile,
					   gp->PolyPTractProfile,
					   allSites->AcceptorSites,
					   l1a,l2a,numAccsites,NUMSITES,Strand,external);

	  sprintf(mess, "U12atac Acceptor Sites \t%8ld", allSites->nAcceptorSites - numAccsites);
	  numAccsites = allSites->nAcceptorSites;
	  printRes(mess);
  }  

  long numDonsites = 0;

  allSites->nDonorSites =
    /* BuildDonors(Sequence,U2,sU2type,sU2, gp->DonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES); */
    BuildDonors(Sequence,U2,sU2type,sU2, gp->DonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
  sprintf (mess,"Donor Sites \t\t%8ld", allSites->nDonorSites);
  numDonsites = allSites->nDonorSites;
  printRes(mess);

  if (U12GTAG){
	  allSites->nDonorSites =
	    BuildDonors(Sequence,U12gtag,sU12type,sU12gtag, gp->U12gtagDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U12gtag Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }
  if (U12ATAC){
	  allSites->nDonorSites =
	    BuildDonors(Sequence, U12atac,sU12type,sU12atac, gp->U12atacDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U12atac Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }
  if (U2GCAG){
	  allSites->nDonorSites =
	    BuildDonors(Sequence,U2, sU2type,sU2gcag, gp->U2gcagDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U2gcag Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTA){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gta, gp->U2gtaDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U2gta Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTG){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gtg, gp->U2gtgDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U2gtg Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  
  if (U2GTY){
	  allSites->nDonorSites =
    	BuildDonors(Sequence,U2,sU2type,sU2gty, gp->U2gtyDonorProfile,allSites->DonorSites,l1b,l2b,numDonsites,NUMSITES,Strand,external);
	  sprintf (mess,"U2gty Donor Sites \t%8ld", allSites->nDonorSites - numDonsites);
	  numDonsites = allSites->nDonorSites;
	  printRes(mess);
  }  

  allSites->nStopCodons =
	GetStopCodons(Sequence,gp->StopProfile, allSites->StopCodons,l1c,l2c);
  sprintf (mess,"Stop Codons \t\t%8ld", allSites->nStopCodons);
  printRes(mess);
  
  if ( U12GTAG || U12ATAC || U2GCAG || U2GTA || U2GTG || U2GTY ){
    /* Predicted sites must be sorted by position */
    printMess ("Sorting Donor and Acceptor sites ...");
    SortSites(allSites->DonorSites,allSites->nDonorSites,donorsites,l1b,l2b);
    SortSites(allSites->AcceptorSites,allSites->nAcceptorSites,acceptorsites,l1a,l2a);
  }
  allSites->nTS=0;
  allSites->nTE=0;
  if (UTR){
    allSites->nTS =
      GetTSS(allSites->TS,allSites->AcceptorSites, allSites->nAcceptorSites, external,hsp,Strand,LengthSequence,l1,l2);
    sprintf(mess, "TS \t\t\t%8ld", allSites->nTS);
    printRes(mess);
    long numTE = 0;
    if(PAS){allSites->nTE =
	GetSitesWithProfile(Sequence,gp->PolyASignalProfile,allSites->TE,l1,l2);
      sprintf(mess, "PolyA Signals \t\t%8ld", allSites->nTE);
      numTE = allSites->nTE;
      printRes(mess);
    }
    allSites->nTE =
      GetTES(allSites->TE,allSites->DonorSites, allSites->nDonorSites,external,hsp,Strand,LengthSequence,l1,l2,numTE);
    sprintf(mess, "TE \t\t\t%8ld", allSites->nTE);
    printRes(mess);
  }

  /* Total number of predicted splice sites in this strand */
  allSites->nSites =
	allSites->nStartCodons +
	allSites->nAcceptorSites +
	allSites->nDonorSites +	
	allSites->nStopCodons + 
        allSites->nTS +
        allSites->nTE;

  sprintf(mess,"---------\t\t%8ld", allSites->nSites);
  printRes(mess);
  

  if ( UTR ){
    /* Predicted sites must be sorted by position */
    printMess ("Sorting TSS/TES sites ...");
    SortSites(allSites->TS,allSites->nTS,tssites,l1,l2);
    SortSites(allSites->TE,allSites->nTE,tesites,l1,l2);
  }
  if (GENAMIC || (!GENAMIC && (EFP || EIP || ETP || ESP || EOP || EXP))){
    /* 2. Building exons with splice sites predicted before */ 
    printMess ("Computing exons ...");   
  

    allExons->nInitialExons =
      BuildInitialExons(allSites->StartCodons,allSites->nStartCodons,
			allSites->DonorSites,allSites->nDonorSites,
			allSites->StopCodons,allSites->nStopCodons,
			gp->MaxDonors,sFIRST,Sequence,
			allExons->InitialExons,NUMEXONS);
    sprintf(mess,"Initial Exons \t\t%8ld", allExons->nInitialExons);
    printRes(mess); 

    allExons->nInternalExons =
      BuildInternalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
			 allSites->DonorSites,allSites->nDonorSites,
			 allSites->StopCodons,allSites->nStopCodons,
			 gp->MaxDonors,sINTERNAL,Sequence,
			 allExons->InternalExons,NUMEXONS);
    sprintf(mess,"Internal Exons \t\t%8ld", allExons->nInternalExons);
    printRes(mess); 

    if (RSS){
      allExons->nZeroLengthExons =
	BuildZeroLengthExons(allSites->AcceptorSites,allSites->nAcceptorSites,
			     allSites->DonorSites,allSites->nDonorSites,
			     allSites->StopCodons,allSites->nStopCodons,
			     gp->MaxDonors,sZEROLENGTH,Sequence,
			     allExons->ZeroLengthExons,NUMEXONS);
      sprintf(mess,"Zero-Length Exons \t%8ld", allExons->nZeroLengthExons);
      printRes(mess); 
    }
    allExons->nTerminalExons =
      BuildTerminalExons(allSites->AcceptorSites,allSites->nAcceptorSites,
			 allSites->StopCodons,allSites->nStopCodons,
			 LengthSequence,cutPoint,sTERMINAL,Sequence,
			 allExons->TerminalExons,NUMEXONS);
    sprintf(mess,"Terminal Exons \t\t%8ld", allExons->nTerminalExons);
    printRes(mess); 
  
    allExons->nSingles =
      BuildSingles(allSites->StartCodons,allSites->nStartCodons,
		   allSites->StopCodons,allSites->nStopCodons,
		   cutPoint, Sequence,
		   allExons->Singles);
    sprintf(mess,"Single genes \t\t%8ld", allExons->nSingles);
    printRes(mess); 
    
    if (UTR){
      allExons->nUtrInitialExons =
	BuildUTRExons(allSites->TS,allSites->nTS,
		      allSites->DonorSites,allSites->nDonorSites,
		      MAXUTRDONORS,MAXUTREXONLENGTH,cutPoint,sUTRFIRST,
		      allExons->UtrInitialExons,NUMEXONS);
      sprintf(mess,"UTR Initial Exons \t%8ld", allExons->nUtrInitialExons);
      printRes(mess);

      allExons->nUtrInitialHalfExons =
	BuildUTRExons(allSites->TS,allSites->nTS,
		      allSites->StartCodons,allSites->nStartCodons,
		      MAXUTRDONORS,MAXUTREXONLENGTH,cutPoint,sUTRFIRSTHALF,
		      allExons->UtrInitialHalfExons,NUMEXONS);
      sprintf(mess,"UTR Initial Half Exons \t%8ld", allExons->nUtrInitialHalfExons);
      printRes(mess);

      allExons->nUtrInternalExons =
	BuildUTRExons(allSites->AcceptorSites,allSites->nAcceptorSites,
		      allSites->DonorSites,allSites->nDonorSites,
		      MAXUTRDONORS,MAXUTREXONLENGTH,cutPoint,sUTRINTERNAL,
		      allExons->UtrInternalExons,NUMEXONS);
      sprintf(mess,"UTR Internal Exons \t%8ld", allExons->nUtrInternalExons);
      printRes(mess); 

      allExons->nUtr5InternalHalfExons =
	BuildUTRExons(allSites->AcceptorSites,allSites->nAcceptorSites,
		      allSites->StartCodons,allSites->nStartCodons,
		      MAXUTRDONORS,MAXUTREXONLENGTH,cutPoint,sUTR5INTERNALHALF,
		      allExons->Utr5InternalHalfExons,NUMEXONS);
      sprintf(mess,"UTR 5' Int. Half Exons \t%8ld", allExons->nUtr5InternalHalfExons);
      printRes(mess);  

      allExons->nUtr3InternalHalfExons =
	BuildUTRExons(allSites->StopCodons,allSites->nStopCodons,
		      allSites->DonorSites,allSites->nDonorSites,
		      MAXUTRDONORS,MAXNMDLENGTH,cutPoint,sUTR3INTERNALHALF,
		      allExons->Utr3InternalHalfExons,NUMEXONS);
      sprintf(mess,"UTR 3' Int. Half Exons \t%8ld", allExons->nUtr3InternalHalfExons);
      printRes(mess);  

      allExons->nUtrTerminalHalfExons =
	BuildUTRExons(allSites->StopCodons,allSites->nStopCodons,
		      allSites->TE,allSites->nTE,
		      MAXUTRDONORS,MAX3UTREXONLENGTH,cutPoint,sUTRTERMINALHALF,
		      allExons->UtrTerminalHalfExons,NUMEXONS);
      sprintf(mess,"UTR Term. Half Exons \t%8ld", allExons->nUtrTerminalHalfExons);
      printRes(mess);   

      allExons->nUtrTerminalExons =
	BuildUTRExons(allSites->AcceptorSites,allSites->nAcceptorSites,
		      allSites->TE,allSites->nTE,
		      MAXUTRDONORS,MAX3UTREXONLENGTH,cutPoint,sUTRTERMINAL,
		      allExons->UtrTerminalExons,NUMEXONS);
      sprintf(mess,"UTR Terminal Exons \t%8ld", allExons->nUtrTerminalExons);
      printRes(mess); 
 
    }

    if (scanORF)
      {
	allExons->nORFs =
	  BuildORFs(allSites->StopCodons,allSites->nStopCodons,
		    allSites->StopCodons,allSites->nStopCodons,
		    cutPoint, Sequence,
		    allExons->ORFs);
	sprintf(mess,"ORFs \t\t\t%8ld", allExons->nORFs);
	printRes(mess); 
      }
    else
      allExons->nORFs = 0;

    /* 3. Scoring and Filtering Exons */
    ScoreExons(Sequence, allExons, 
	       l1, l2, Strand, 
	       external, hsp,
	       isochores,nIsochores,
	       GCInfo);
  
    /* Total number of built exons in this strand */
    allExons->nExons =
      allExons->nInitialExons +
      allExons->nInternalExons +
      allExons->nZeroLengthExons +
      allExons->nTerminalExons +
      allExons->nSingles +
      allExons->nORFs +
      allExons->nUtrInitialExons +
      allExons->nUtrInitialHalfExons +
      allExons->nUtrInternalExons +
      allExons->nUtr5InternalHalfExons +
      allExons->nUtr3InternalHalfExons +
      allExons->nUtrTerminalExons +
      allExons->nUtrTerminalHalfExons;

    sprintf(mess,"---------\t\t%8ld", allExons->nExons);
    printRes(mess); 
  }
}
