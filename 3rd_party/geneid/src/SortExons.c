/*************************************************************************
*                                                                        *
*   Module: SortExons                                                    *
*                                                                        *
*   Sort by left signal position all of predicted exons                  *
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

/*  $Id: SortExons.c,v 1.18 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum allowed number of generic exons (multiplied by FSORT) */
extern long NUMEXONS;
extern int RSS;
extern int EVD;
extern int UTR;
extern int FWD,RVS;
extern int scanORF;

/* Artificial initial gene feature: force complete gene prediction */
void InsertBeginExon(exonGFF* Exons, long lowerlimit)
{
  exonGFF* e;
  exonGFF* re;

  /* 1. Forward Strand */
  /* Create the exon structure */
  if ((e = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: artificial FWD Begin exon");

  /* Create the sites structure */
  if ((e->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for artificial FWD Begin exon");

  /* Create the sites structure */
  if ((e->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for artificial FWD Begin exon");
  
  e->Acceptor->Position = lowerlimit;
  e->Donor->Position = lowerlimit;
  e->Acceptor->class = U2;
  strcpy(e->Acceptor->subtype,sU2);
  strcpy(e->Acceptor->type,sU2type);
  e->Donor->class = U2;
  strcpy(e->Donor->subtype,sU2);
  strcpy(e->Donor->type,sU2type);

  /* Save the extracted exon */
  Exons[0].Acceptor = e->Acceptor;
  Exons[0].Donor = e->Donor;
  strcpy(Exons[0].Type,sBEGIN);
  Exons[0].Frame = 0;
  Exons[0].Strand = '+';
  Exons[0].Score = MAXSCORE;
  Exons[0].PartialScore = 0.0;
  Exons[0].HSPScore = 0.0;
  Exons[0].R = 0.0;
  Exons[0].GeneScore = 0.0; 
  Exons[0].Remainder = 0;
  strcpy(Exons[0].Group,NOGROUP);
  Exons[0].evidence = 0;
  Exons[0].offset1 = 0;
  Exons[0].offset2 = 0; 
  Exons[0].lValue = 0;
  Exons[0].rValue = 0; 
  Exons[0].selected = 0;

  /* 2. Reverse Strand */
  /* Create the exon structure */
  if ((re = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: artificial RVS Begin exon");

  /* Create the sites structure */
  if ((re->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for artificial RVS Begin exon");
  /* Create the sites structure */
  if ((re->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for artificial RVS Begin exon");
  
  re->Acceptor->Position = 0;
  re->Donor->Position = 0;
  re->Acceptor->class = U2;
  strcpy(re->Acceptor->subtype,sU2);
  strcpy(re->Acceptor->type,sU2type);
  re->Donor->class = U2;
  strcpy(re->Donor->subtype,sU2);
  strcpy(re->Donor->type,sU2type);

  /* Save the extracted exon */
  Exons[1].Acceptor = re->Acceptor;
  Exons[1].Donor = re->Donor;
  strcpy(Exons[1].Type,sBEGIN);
  Exons[1].Frame = 0;
  Exons[1].Strand = '-';
  Exons[1].Score = MAXSCORE;
  Exons[1].PartialScore = 0.0;
  Exons[1].HSPScore = 0.0;
  Exons[1].R = 0.0;
  Exons[1].GeneScore = 0.0; 
  Exons[1].Remainder = 0;
  strcpy(Exons[1].Group,NOGROUP);
  Exons[1].evidence = 0;
  Exons[1].offset1 = 0;
  Exons[1].offset2 = 0; 
  Exons[1].lValue = 0;
  Exons[1].rValue = 0; 
  Exons[1].selected = 0;

}

/* Artificial terminal gene feature: force complete gene prediction */
void InsertEndExon(exonGFF* Exons, long n, long L)
{
  exonGFF* e;
  exonGFF* re;

  /* 1. Forward Strand */
  /* Create the exon structure */
  if ((e = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: artificial FWD End exon");

  /* Create the sites structure */
  if ((e->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for artificial FWD End exon");

  /* Create the sites structure */
  if ((e->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for artificial FWD End exon");

  e->Acceptor->Position = L;
  e->Donor->Position = L;
  e->Acceptor->class = U2;
  strcpy(e->Acceptor->subtype,sU2);
  strcpy(e->Acceptor->type,sU2type);
  e->Donor->class = U2;
  strcpy(e->Donor->subtype,sU2);
  strcpy(e->Donor->type,sU2type);

  /* Save the extracted exon */
  Exons[n].Acceptor = e->Acceptor;
  Exons[n].Donor = e->Donor;
  strcpy(Exons[n].Type,sEND);
  Exons[n].Frame = 0;
  Exons[n].Strand = '+';
  Exons[n].Score = MAXSCORE;
  Exons[n].PartialScore = 0.0;
  Exons[n].HSPScore = 0.0;
  Exons[n].R = 0.0;
  Exons[n].GeneScore = 0.0; 
  Exons[n].Remainder = 0;
  strcpy(Exons[n].Group,NOGROUP);
  Exons[n].evidence = 0;
  Exons[n].offset1 = 0;
  Exons[n].offset2 = 0; 
  Exons[n].lValue = 0;
  Exons[n].rValue = 0; 
  Exons[n].selected = 0;

  /* 1. Reverse Strand */
  /* Create the exon structure */
  if ((re = (exonGFF*) malloc(sizeof(exonGFF))) == NULL)
        printError("Not enough memory: artificial RVS End exon");

  /* Create the sites structure */
  if ((re->Acceptor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor site for artificial RVS End exon");

  /* Create the sites structure */
  if ((re->Donor =
       (struct s_site *) malloc(sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor site for artificial RVS End exon");

  re->Acceptor->Position = L;
  re->Donor->Position = L;
  re->Acceptor->class = U2;
  strcpy(re->Acceptor->subtype,sU2);
  strcpy(re->Acceptor->type,sU2type);
  re->Donor->class = U2;
  strcpy(re->Donor->subtype,sU2);
  strcpy(re->Donor->type,sU2type);

  /* Save the extracted exon */
  Exons[n+1].Acceptor = re->Acceptor;
  Exons[n+1].Donor = re->Donor;
  strcpy(Exons[n+1].Type,sEND);
  Exons[n+1].Frame = 0;
  Exons[n+1].Strand = '-';
  Exons[n+1].Score = MAXSCORE;
  Exons[n+1].PartialScore = 0.0;
  Exons[n+1].HSPScore = 0.0;
  Exons[n+1].R = 0.0;
  Exons[n+1].GeneScore = 0.0; 
  Exons[n+1].Remainder = 0;
  strcpy(Exons[n+1].Group,NOGROUP);
  Exons[n+1].evidence = 0;
  Exons[n+1].offset1 = 0;
  Exons[n+1].offset2 = 0; 
  Exons[n+1].lValue = 0;
  Exons[n+1].rValue = 0; 
  Exons[n+1].selected = 0;

}

/* Struct for a node (list): pointer to exon and to next node */
struct exonitem 
{
  exonGFF* Exon;
  struct exonitem* nexitem;
} ExonItem;

/* Set free a list of nodes */
void FreeItems(struct exonitem *q)
{
  if (q == NULL)
    ;
  else
    {
      FreeItems(q->nexitem);
      free(q);
    }
}

/* Insert an exon in the selected list according to its left position */
void UpdateList(struct exonitem** p, exonGFF* InputExon)
{
  /* Insert the new node at the end of the list */
  while (*p!=NULL)
    p=&((*p)->nexitem); 
  
  /* Allocating new node for this exon */
  if ((*p= (struct exonitem *) malloc(sizeof(struct exonitem))) == NULL)
    printError("Not enough memory: new exonitem (sorting)");
  
  /* Updating information for the node */
  (*p)->Exon=InputExon;
  (*p)->nexitem=NULL;
}

/* Sort all of predicted exons by placing them in an array of lists */
/* corresponding every list to a beginning position for predicted exons */
void SortExons(packExons* allExons,
               packExons* allExons_r, 
               packExternalInformation* external,
	       packEvidence* pv,
	       exonGFF* Exons,         
               long l1, long l2,long lowerlimit,
	       long upperlimit)
{ 
  struct exonitem **ExonList, *q;
  long i;
  long acceptor;
  /* long donor; */
  long n;
  int offset;
  long l;
  long left;
  long right;
  long room;
  char mess[MAXSTRING];
  long HowMany;

  /* 0. Creating the array for sorting: 1 - Length of fragment */
  left = l1 + COFFSET - LENGTHCODON;
  right = l2 + COFFSET + LENGTHCODON;
  l =  right - left + 1;

  /* Allocating memory for array ExonList */
  if ((ExonList = 
       (struct exonitem **) calloc(l, sizeof(struct exonitem *))) == NULL)
    printError("Not enough memory: ExonList array (sorting)");

  /* Reset the positions, pointing to NULL */
  for (i=0;i<l;i++)
    ExonList[i]=NULL;

  /* 1. Insert exons in the proper list according to its left position */
  /* Adding predicted exons in the forward sense */
  for (i=0; i<allExons->nInitialExons; i++) 
    {
      /* Correction between real and relative to fragment position */
      acceptor=(allExons->InitialExons+i)->Acceptor->Position - left;
      (allExons->InitialExons+i)->Strand = '+';
      /* Fixing positions according to the COFFSET in arrays */
      CorrectExon(allExons->InitialExons+i);
      offset = (allExons->InitialExons+i)->offset1;
      /* Insert exon in the proper list depending on the left position */
      UpdateList(&(ExonList[acceptor+offset]), allExons->InitialExons+i);
    }

  for (i=0;i<allExons->nInternalExons;i++) 
    {
      acceptor=(allExons->InternalExons+i)->Acceptor->Position - left;
      (allExons->InternalExons+i)->Strand = '+';
      CorrectExon(allExons->InternalExons+i);
      offset = (allExons->InternalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->InternalExons+i); 
    }
  if (RSS){
    for (i=0;i<allExons->nZeroLengthExons;i++) 
      {
	acceptor=(allExons->ZeroLengthExons+i)->Acceptor->Position - left;
	(allExons->ZeroLengthExons+i)->Strand = '+';
	CorrectExon(allExons->ZeroLengthExons+i);
	offset = (allExons->ZeroLengthExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->ZeroLengthExons+i); 
      }
  }
  for (i=0;i<allExons->nTerminalExons;i++) 
    {
      acceptor=(allExons->TerminalExons+i)->Acceptor->Position - left;
      (allExons->TerminalExons+i)->Strand = '+';
      CorrectExon(allExons->TerminalExons+i);
      offset = (allExons->TerminalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->TerminalExons+i); 
    }
    
  for (i=0;i<allExons->nSingles;i++) 
    {
      acceptor=(allExons->Singles+i)->Acceptor->Position - left;
      (allExons->Singles+i)->Strand = '+';
      CorrectExon(allExons->Singles+i);
      offset = (allExons->Singles+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons->Singles+i); 
    }
  
  if (scanORF)
    for (i=0;i<allExons->nORFs;i++) 
      {
	acceptor=(allExons->ORFs+i)->Acceptor->Position - left;
	(allExons->ORFs+i)->Strand = '+';
	CorrectORF(allExons->ORFs+i);
	offset = (allExons->ORFs+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->ORFs+i); 
      }
  if (UTR){
    for (i=0;i<allExons->nUtrInitialExons;i++) 
      {
	acceptor=(allExons->UtrInitialExons+i)->Acceptor->Position - left;
	(allExons->UtrInitialExons+i)->Strand = '+';
	CorrectExon(allExons->UtrInitialExons+i);
	offset = (allExons->UtrInitialExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->UtrInitialExons+i); 
      }
    for (i=0;i<allExons->nUtrInitialHalfExons;i++)
      {
	acceptor=(allExons->UtrInitialHalfExons+i)->Acceptor->Position - left;
	(allExons->UtrInitialHalfExons+i)->Strand = '+';
	CorrectUTR(allExons->UtrInitialHalfExons+i);
	offset = (allExons->UtrInitialHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->UtrInitialHalfExons+i);
      }
    for (i=0;i<allExons->nUtrInternalExons;i++)
      {
	acceptor=(allExons->UtrInternalExons+i)->Acceptor->Position - left;
	(allExons->UtrInternalExons+i)->Strand = '+';
	CorrectExon(allExons->UtrInternalExons+i);
	offset = (allExons->UtrInternalExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->UtrInternalExons+i);
      }
    for (i=0;i<allExons->nUtr5InternalHalfExons;i++)
      {
	acceptor=(allExons->Utr5InternalHalfExons+i)->Acceptor->Position - left;
	(allExons->Utr5InternalHalfExons+i)->Strand = '+';
	CorrectUTR(allExons->Utr5InternalHalfExons+i);
	offset = (allExons->Utr5InternalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->Utr5InternalHalfExons+i);
      }
    for (i=0;i<allExons->nUtr3InternalHalfExons;i++)
      {
	acceptor=(allExons->Utr3InternalHalfExons+i)->Acceptor->Position - left;
	(allExons->Utr3InternalHalfExons+i)->Strand = '+';
	CorrectUTR(allExons->Utr3InternalHalfExons+i);
	offset = (allExons->Utr3InternalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->Utr3InternalHalfExons+i);
      }
    for (i=0;i<allExons->nUtrTerminalHalfExons;i++)
      {
	acceptor=(allExons->UtrTerminalHalfExons+i)->Acceptor->Position - left;
	(allExons->UtrTerminalHalfExons+i)->Strand = '+';
	CorrectUTR(allExons->UtrTerminalHalfExons+i);
	offset = (allExons->UtrTerminalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->UtrTerminalHalfExons+i);
      }
    for (i=0;i<allExons->nUtrTerminalExons;i++)
      {
	acceptor=(allExons->UtrTerminalExons+i)->Acceptor->Position - left;
	(allExons->UtrTerminalExons+i)->Strand = '+';
	CorrectExon(allExons->UtrTerminalExons+i);
	offset = (allExons->UtrTerminalExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons->UtrTerminalExons+i);
      }
  }
  /* Adding predicted exons in the reverse sense */
  for (i=0;i<allExons_r->nInitialExons;i++) 
    {
      acceptor=(allExons_r->InitialExons+i)->Acceptor->Position - left;
      (allExons_r->InitialExons+i)->Strand = '-';
      CorrectExon(allExons_r->InitialExons+i);
      offset = (allExons_r->InitialExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->InitialExons+i); 
    }
  for (i=0;i<allExons_r->nInternalExons;i++) 
    {
      acceptor=(allExons_r->InternalExons+i)->Acceptor->Position - left;
      (allExons_r->InternalExons+i)->Strand = '-';
      CorrectExon(allExons_r->InternalExons+i);
      offset = (allExons_r->InternalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->InternalExons+i); 
    }

  if (RSS){
    for (i=0;i<allExons_r->nZeroLengthExons;i++) 
      {
	acceptor=(allExons_r->ZeroLengthExons+i)->Acceptor->Position - left;
	(allExons_r->ZeroLengthExons+i)->Strand = '-';
	CorrectExon(allExons_r->ZeroLengthExons+i);
	offset = (allExons_r->ZeroLengthExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->ZeroLengthExons+i); 
      }
  }

  for (i=0;i<allExons_r->nTerminalExons;i++) 
    {   
      acceptor=(allExons_r->TerminalExons+i)->Acceptor->Position - left;
      (allExons_r->TerminalExons+i)->Strand = '-';
      CorrectExon(allExons_r->TerminalExons+i);
      offset = (allExons_r->TerminalExons+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->TerminalExons+i); 
    }

  for (i=0;i<allExons_r->nSingles;i++) 
    {
      acceptor=(allExons_r->Singles+i)->Acceptor->Position - left;
      (allExons_r->Singles+i)->Strand = '-';
      CorrectExon(allExons_r->Singles+i);
      offset = (allExons_r->Singles+i)->offset1;
      UpdateList(&(ExonList[acceptor + offset]), allExons_r->Singles+i); 
    }

  if (scanORF)
    for (i=0;i<allExons_r->nORFs;i++) 
      {
	acceptor=(allExons_r->ORFs+i)->Acceptor->Position - left;
	(allExons_r->ORFs+i)->Strand = '-';
	CorrectORF(allExons_r->ORFs+i);
	offset = (allExons_r->ORFs+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->ORFs+i); 
      }

  if (UTR){
    for (i=0;i<allExons_r->nUtrInitialExons;i++) 
      {
	acceptor=(allExons_r->UtrInitialExons+i)->Acceptor->Position - left;
	(allExons_r->UtrInitialExons+i)->Strand = '-';
	CorrectExon(allExons_r->UtrInitialExons+i);
	offset = (allExons_r->UtrInitialExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->UtrInitialExons+i); 
      }
    for (i=0;i<allExons_r->nUtrInitialHalfExons;i++)
      {
	acceptor=(allExons_r->UtrInitialHalfExons+i)->Acceptor->Position - left;
	(allExons_r->UtrInitialHalfExons+i)->Strand = '-';
	CorrectUTR(allExons_r->UtrInitialHalfExons+i);
	offset = (allExons_r->UtrInitialHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->UtrInitialHalfExons+i);
      }
    for (i=0;i<allExons_r->nUtrInternalExons;i++)
      {
	acceptor=(allExons_r->UtrInternalExons+i)->Acceptor->Position - left;
	(allExons_r->UtrInternalExons+i)->Strand = '-';
	CorrectExon(allExons_r->UtrInternalExons+i);
	offset = (allExons_r->UtrInternalExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->UtrInternalExons+i);
      }
    for (i=0;i<allExons_r->nUtr5InternalHalfExons;i++)
      {
	acceptor=(allExons_r->Utr5InternalHalfExons+i)->Acceptor->Position - left;
	(allExons_r->Utr5InternalHalfExons+i)->Strand = '-';
	CorrectUTR(allExons_r->Utr5InternalHalfExons+i);
	offset = (allExons_r->Utr5InternalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->Utr5InternalHalfExons+i);
      }
    for (i=0;i<allExons_r->nUtr3InternalHalfExons;i++)
      {
	acceptor=(allExons_r->Utr3InternalHalfExons+i)->Acceptor->Position - left;
	(allExons_r->Utr3InternalHalfExons+i)->Strand = '-';
	CorrectUTR(allExons_r->Utr3InternalHalfExons+i);
	offset = (allExons_r->Utr3InternalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->Utr3InternalHalfExons+i);
      }
    for (i=0;i<allExons_r->nUtrTerminalHalfExons;i++)
      {
	acceptor=(allExons_r->UtrTerminalHalfExons+i)->Acceptor->Position - left;
	(allExons_r->UtrTerminalHalfExons+i)->Strand = '-';
	CorrectUTR(allExons_r->UtrTerminalHalfExons+i);
	offset = (allExons_r->UtrTerminalHalfExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->UtrTerminalHalfExons+i);
      }
    for (i=0;i<allExons_r->nUtrTerminalExons;i++)
      {
	acceptor=(allExons_r->UtrTerminalExons+i)->Acceptor->Position - left;
	(allExons_r->UtrTerminalExons+i)->Strand = '-';
	CorrectExon(allExons_r->UtrTerminalExons+i);
	offset = (allExons_r->UtrTerminalExons+i)->offset1;
	UpdateList(&(ExonList[acceptor + offset]), allExons_r->UtrTerminalExons+i);
      }
  }
  /* 2. Merge evidence exons (annotations) with predictions */
  if (EVD && pv != NULL)
    for (i = external->i1vExons; i < external->i2vExons ; i++) 
      {  
	acceptor=(pv->vExons + i)->Acceptor->Position - left;
	offset = (pv->vExons + i)->offset1;
	room = ((acceptor + offset)<0)? 0 : (acceptor + offset);
		
	/* Requirement 1: range of values */
	if (acceptor < 0)
	  {
	    /* Skip wrong exon (A) */
	    sprintf(mess,"Skipped evidence (range): %ld %ld %c %hd",
		    (pv->vExons + i)->Acceptor->Position,
		    (pv->vExons + i)->Donor->Position,
		    (pv->vExons + i)->Strand,
		    (pv->vExons + i)->Frame);
	    printMess(mess);
	  }
	else
	  {
	    /* Requirement 2: forced strand prediction */
	    if ( (!FWD && (pv->vExons + i)->Strand == '+') ||
		 (!RVS && (pv->vExons + i)->Strand == '-'))
	      {
		/* Skip wrong exon (B) */
		sprintf(mess,"Skipped evidence (strand): %ld %ld %c %hd",
			(pv->vExons + i)->Acceptor->Position,
			(pv->vExons + i)->Donor->Position,
			(pv->vExons + i)->Strand,
			(pv->vExons + i)->Frame);
		printMess(mess);
	      }
	    else{
	      UpdateList(&(ExonList[room]), pv->vExons+i);
	    }
	  }
      }

  /* FIRST SPLIT: Insert first artificial exon */
  if (l1 == lowerlimit)
    {
      InsertBeginExon(Exons,lowerlimit);
      n = 2;
    }
  else
    n = 0;

  /* 3. Traversing the table extracting the exons sorted by left position */
  HowMany = FSORT*NUMEXONS;
  /* for (i=0, n=0; i<l; i++) */
  for (i=0; i<l; i++)
    {
      /* Extracting exons from list q */
      q=ExonList[i];
      while (q != NULL)
		{
		  /* Save the extracted exon */
		  Exons[n].Acceptor = q->Exon->Acceptor;
		  Exons[n].Donor = q->Exon->Donor;
		  strcpy(Exons[n].Type,q->Exon->Type);
		  Exons[n].Frame = q->Exon->Frame;
		  Exons[n].Strand = q->Exon->Strand;
		  Exons[n].Score = q->Exon->Score;
		  Exons[n].PartialScore = q->Exon->PartialScore;
		  Exons[n].HSPScore = q->Exon->HSPScore;
		  Exons[n].R = q->Exon->R;
		  Exons[n].GeneScore = 0.0; 
		  Exons[n].Remainder = q->Exon->Remainder;
		  strcpy(Exons[n].Group,q->Exon->Group);
		  Exons[n].evidence = q->Exon->evidence;
		  Exons[n].offset1 = q->Exon->offset1;
		  Exons[n].offset2 = q->Exon->offset2; 
		  Exons[n].lValue = q->Exon->lValue;
		  Exons[n].rValue = q->Exon->rValue; 
		  Exons[n].selected = 0;

		  n++;
	  
		  if (n >= HowMany)
			printError("Too many predicted exons: increase FSORT parameter");
	  
		  q=q->nexitem;
		}

	  /* Free chained items in the processed list q */
	  FreeItems(ExonList[i]);       
	}

  /* FINISHING split: Insert last artificial exon */
  if (l2 == upperlimit)
	{
	  InsertEndExon(Exons, n, upperlimit + 1);
	  n = n + 2;
	  
	  if (n >= HowMany)
		printError("Too many predicted exons: increase FSORT parameter");
	}
  
  /* Free empty array */
  free(ExonList);
}

   
