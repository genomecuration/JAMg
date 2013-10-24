/*************************************************************************
*                                                                        *
*   Module: Translate                                                    *
*                                                                        *
*   Translate genomic sequences into protein products (amino acids)      *
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

/*  $Id: Translate.c,v 1.15 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Obtaining the amino acid sequence from an exon (given a reading frame) */
int Translate(long p1,
              long p2,
              short fra,
              short rmd,
              char* s,
              dict* dAA,
              char sAux[])
{
  int nAA, nAux;
  long i;
  char codon[LENGTHCODON+1];

  /* Saving first uncomplete codon in the exon (frame) */
  for(i=0, nAux = 0; i<fra; i++)
	{
	  sAux[i] = s[p1+i] + 32;
	  nAux++;
	}

  /* Translating complet codons in the exon */
  nAA = nAux;
  for(i=p1+ fra; i<= p2 - rmd; i=i+3)
	{
	  codon[0] = s[i];
	  codon[1] = s[i+1];
	  codon[2] = s[i+2];
	  codon[3] = '\0';
        
	  /* Looking up the amino acid dictionary */
	  sAux[nAux] = getAADict(dAA,codon);
	  nAux++;
	}
  nAA = nAux - nAA;
   
  /* Saving last uncomplete codon in the exon (remainder) */
  for(i=p2 - rmd + 1; i <= p2; i++)
	{
	  sAux[nAux] = s[i] + 32;
	  nAux++;
	}
  sAux[nAux] = '\0';

  /* Adding to the count both uncomplete codons */
  if (fra)
	nAA++;
  if (rmd)
	nAA++;

  return(nAA);
}

/* Translate gene sequence into the protein */
void TranslateGene(exonGFF* e,
                   char* s,
                   dict* dAA,
                   long nExons,
                   int** tAA,
                   char* prot,
                   long* nAA)
{  
  long i;
  short j;
  int totalAA;
  int currAA;
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  char codon[LENGTHCODON+1];
  char aa;
  short currFrame;
  short currRmd;
  int lAux;
  char rmdProt[LENGTHCODON+1];
  int lastExon = 1;

  /* A. Translating a forward sense exon: Terminal > Internal >.. First */
  if (e->Strand == '+')
    {
      prot[0] = '\0';
      /* Traversing the list of exons */      
      for(i=0, totalAA=0, currFrame=0; i<nExons; i++)
	{
	  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){
	    /* 0. Acquire the real positions of the exon in the sequence */
	    p1 = (e->evidence)?
	      e->Acceptor->Position + e->offset1 : 
	      e->Acceptor->Position + e->offset1 - COFFSET;
	    p2 = (e->evidence)? 
	      e->Donor->Position + e->offset2: 
	      e->Donor->Position + e->offset2 - COFFSET;

	    /* -- Remainder nucleotides of the last exon in the gene -- */
	    if (!i || lastExon)
	      {
		switch (e->Remainder)
		  {
		  case 0:
		    rmdProt[0] = '\0'; 
		    break;
		  case 1:
		    rmdProt[0] = s[p2] + 32; 
		    rmdProt[1] = '\0'; 
		    break;
		  case 2:
		    /* curr RMD must be TWO */
		    rmdProt[0] = s[p2-1] + 32;
		    rmdProt[1] = s[p2] + 32;
		    rmdProt[2] = '\0'; 
		    break;
		  } 
		/* Leaving the result into prot */
		sprintf(prot,"%s",rmdProt);
	      }
	  
	    /* 1. Translating current exon sequence */
	    currAA = Translate(p1 + e->Frame,
			       p2 - (3 - e->Remainder)%3,
			       0,
			       0,
			       s, dAA, sAux);

	    /* 2. Translating shared codon between current and previous exon */
	    /* concat after the translation of the current exon in sAux */
	    /* this amino acid has been counted into Translate */
	    switch (currFrame)
	      {
	      case 1:
		/* curr RMD must be TWO */
		codon[0] = s[p2-1];
		codon[1] = s[p2];
		codon[3] = '\0';
		/* translate this codon */
		aa = getAADict(dAA,codon);
		lAux = strlen(sAux);
		sAux[lAux] = aa;
		sAux[lAux+1] = '\0'; 
		break;
            
	      case 2:
		/* curr RMD must be ONE */
		codon[0] = s[p2];
		codon[3] = '\0';
		aa = getAADict(dAA,codon);
		lAux = strlen(sAux);
		sAux[lAux] = aa;
		sAux[lAux+1] = '\0'; 
		break;
	      }

	    /* Concat translation(current exon) with translation(last exon) */
	    strcat(sAux,prot);
	    /* Leaving the result into prot */
	    strcpy(prot,sAux);

	    /* 3. Saving information to translate the next shared codon */
	    currFrame = e->Frame;
	    switch (currFrame)
	      {
	      case 1:
		codon[2] = s[p1];
		break;
            
	      case 2:
		/* curr RMD must be TWO */
		codon[1] = s[p1];
		codon[2] = s[p1+1];
		break;
	      }

	    /* 4. Updating amino acids range for this exon. begin(0), end(1) */
	    tAA[i][0] = MAX(1,totalAA);
        
	    /* Discounting one uncomplete codon: twice added */
	    if (!e->Remainder && i && !lastExon)
	      tAA[i][0]++;
           
	    if (e->Frame)
	      totalAA++;
        
	    totalAA = totalAA + currAA;
	    tAA[i][1] = totalAA;
	  }
	  /* 5. Pointer jumping to the next exon */
	  e = e->PreviousExon;
	  lastExon = 0;
	} /* endfor */
      
	  /* 6. Adding first uncomplete codon of the first exon in the gene */    
      switch (currFrame)
	{
	case 0:
	  rmdProt[0] = '\0'; 
	  break;
	case 1:
	  rmdProt[0] = codon[2] + 32; 
	  rmdProt[1] = '\0'; 
	  break;
	case 2:
	  /* curr RMD must be TWO */
	  rmdProt[0] = codon[1] + 32;
	  rmdProt[1] = codon[2] + 32;
	  rmdProt[2] = '\0'; 
	  break;
	}
      
      /* Concat uncomplete codon at the beginning of the protein */
      sprintf(sAux,"%s%s",rmdProt,prot);
      strcpy(prot,sAux);
		
    }

  /* B. Translating a reverse sense exon: First > Internal >.. Terminal */
  else
    {
      prot[0] = '\0'; 
      /* Traversing the list of exons */
      for(i=0, totalAA=0, currRmd=0; i<nExons; i++)
	{
	  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){
	    p1 = (e->evidence)? 
	      e->Acceptor->Position + e->offset1: 
	      e->Acceptor->Position + e->offset1 - COFFSET;
	    p2 = (e->evidence)? 
	      e->Donor->Position + e->offset2: 
	      e->Donor->Position + e->offset2 - COFFSET;
	    /* Memory for the reverse exon sequence */
	    if ((rs = (char*) calloc(p2-p1+2,sizeof(char))) == NULL)
	      printError("Not enough memory: reverse gene translation");
	 
	    ReverseSubSequence(p1, p2, s, rs);

	    /* Frame of the last exon in sequence */
	    if (!i || lastExon)
	      {
		switch (e->Frame)
		  {
		  case 0:
		    rmdProt[0] = '\0'; 
		    break;
		  case 1:
		    rmdProt[0] = rs[0] + 32; 
		    rmdProt[1] = '\0'; 
		    break;
		  case 2:
		    /* curr RMD must be TWO */
		    rmdProt[0] = rs[0] + 32;
		    rmdProt[1] = rs[1] + 32;
		    rmdProt[2] = '\0'; 
		    break;
		  }	
	      }
	  
	    /* Translating the current exon (forward reading) */
	    currAA = Translate(0 + e->Frame, p2-p1-(3 - e->Remainder)%3,
			       0,
			       0,
			       rs, dAA, sAux);
	  
	    /* Translating shared codon between current and previous exon */
	    switch (currRmd)
	      {
	      case 1:
		/* curr FRAME must be TWO */
		codon[1] = rs[0];
		codon[2] = rs[1];
		codon[3] = '\0';
		/* translate this codon */
		aa = getAADict(dAA,codon);
		lAux = strlen(prot);
		prot[lAux] = aa;
		prot[lAux+1] = '\0'; 
		break;
	      
	      case 2:
		/* curr FRAME must be ONE */
		codon[2] = rs[0];
		codon[3] = '\0';
		aa = getAADict(dAA,codon);
		lAux = strlen(prot);
		prot[lAux] = aa;
		prot[lAux+1] = '\0'; 
		break;
	      }
	  
	    /* exon translation after translated remainder */
	    strcat(prot,sAux);
	  
	    /* Codon information for translating the next shared codon */
	    currRmd = (3 - e->Remainder)%3;
	    switch (currRmd)
	      {
	      case 1:
		codon[0] = rs[p2-p1];
		break;
	      
	      case 2:
		/* curr RMD must be TWO */
		codon[0] = rs[p2-p1-1];
		codon[1] = rs[p2-p1];
		break;
	      }
	  
	    /* Updating boundary aminoacids of this exon */
	    tAA[i][0] = MAX(1,totalAA);

	    if (!e->Frame && i && !lastExon)
	      tAA[i][0]++;
  	  
	    if (currRmd)
	      totalAA++;
	  
	    totalAA = totalAA + currAA;
	    tAA[i][1] = totalAA;
	  
	    /* Next exon */
	    free(rs);
	  }
	  e = e->PreviousExon;
	  lastExon = 0;
		  
	}
      /* Frame of the last exon in the gene */
      sprintf(sAux,"%s%s",rmdProt,prot);
      strcpy(prot,sAux);
      
      /* Remainder of the first exon in sequence */
      switch (currRmd)
	{
	case 0:
	  rmdProt[0] = '\0'; 
	  break;
	case 1:
	  rmdProt[0] = codon[0] + 32; 
	  rmdProt[1] = '\0'; 
	  break;
	case 2:
	  /* curr RMD must be TWO */
	  rmdProt[0] = codon[0] + 32;
	  rmdProt[1] = codon[1] + 32;
	  rmdProt[2] = '\0'; 
	  break;
	}
      
      lAux = strlen(prot);
      for(j = 0; j < currRmd; j++)
	prot[lAux+j] = rmdProt[j];
      prot[lAux+j] = '\0'; 
      
    } /* end of reverse translation */

  *nAA = totalAA;
}

/* Extract the genomic (CDS) sequence of a predicted gene */
/* Returns the length of the genomic sequence produced */
void GetcDNA(exonGFF* e,
             char* s,
             long nExons,
             char* cDNA,
             long* nNN)
{  
  char* tmpDNA;
  long p1,p2;
  char* rs;
  int i;
  long j;

  if ((tmpDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
    printError("Not enough memory: producing CDS DNA");

  cDNA[0] = '\0';

  if (e->Strand == '+')
    {
      for(i=0, *nNN = 0; i<nExons; i++)
		{
		  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){
		    p1 = (e->evidence)? 
		      e->Acceptor->Position + e->offset1: 
		      e->Acceptor->Position + e->offset1 - COFFSET;
		    p2 = (e->evidence)? 
		      e->Donor->Position + e->offset2: 
		      e->Donor->Position + e->offset2 - COFFSET;
	  
		  /* Get the cDNA for this exon */
		  tmpDNA[0] = '\0';
   
		  for(j=p1; j <= p2; j++)
			tmpDNA[j-p1] = s[j];
        
		  tmpDNA[j-p1]='\0';
        
		  /* Concat the current cDNA before the previous */
		  strcat(tmpDNA,cDNA);
		  strcpy(cDNA,tmpDNA);
        
		  *nNN += p2-p1+1;
		  }
		  e = e->PreviousExon;
		}
    }
  /* Reverse and complement the sequence in this strand */
  else
    {
      if ((rs = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
		printError("Not enough memory: producing reverse CDS DNA");

      for(i=0, *nNN = 0; i<nExons; i++)
		{
		  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){
		    p1 = (e->evidence)? 
		      e->Acceptor->Position + e->offset1: 
		      e->Acceptor->Position + e->offset1 - COFFSET;
		    p2 = (e->evidence)? 
		      e->Donor->Position + e->offset2: 
		      e->Donor->Position + e->offset2 - COFFSET;
		  /* Get the cDNA for this exon */
		  tmpDNA[0] = '\0';
		  rs[0] = '\0';
   
		  for(j=p1; j <= p2; j++)
			tmpDNA[j-p1] = s[j];
        
		  tmpDNA[j-p1]='\0';
        
		  ReverseSubSequence(0, j-p1-1, tmpDNA, rs);
		  rs[j-p1]='\0';
   
		  /* Concat the current cDNA after the previous */
		  strcat(cDNA,rs);
   
		  *nNN += p2-p1+1;
		  }
		  e = e->PreviousExon;
		  
		} 
      /* Set free reverse sequence */
      free(rs);
    }

  /* Set free temporary result sequence */
  free(tmpDNA);
}

/* Extract the genomic (exonic) sequence of a predicted gene */
/* Returns the length of the genomic sequence produced */
void GetTDNA(exonGFF* e,
             char* s,
             long nExons,
             char* cDNA,
             long* nNN)
{  
  char* tmpDNA;
  long p1,p2;
  char* rs;
  int i;
  long j;

  if ((tmpDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
    printError("Not enough memory: producing exonic DNA");

  cDNA[0] = '\0';

  if (e->Strand == '+')
    {
      for(i=0, *nNN = 0; i<nExons; i++)
	{
/* 	  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){ */
	    p1 = (e->evidence)? 
	      e->Acceptor->Position + e->offset1: 
	      e->Acceptor->Position + e->offset1 - COFFSET;
	    p2 = (e->evidence)? 
	      e->Donor->Position + e->offset2: 
	      e->Donor->Position + e->offset2 - COFFSET;
	  
	    /* Get the cDNA for this exon */
	    tmpDNA[0] = '\0';
   
	    for(j=p1; j <= p2; j++)
	      tmpDNA[j-p1] = s[j];
        
	    tmpDNA[j-p1]='\0';
        
	    /* Concat the current cDNA before the previous */
	    strcat(tmpDNA,cDNA);
	    strcpy(cDNA,tmpDNA);
        
	    *nNN += p2-p1+1;
	  /* } */
	  e = e->PreviousExon;
	}
    }
  /* Reverse and complement the sequence in this strand */
  else
    {
      if ((rs = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
	printError("Not enough memory: producing reverse exonic DNA");

      for(i=0, *nNN = 0; i<nExons; i++)
	{
/* 	  if (!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)){ */
	    p1 = (e->evidence)? 
	      e->Acceptor->Position + e->offset1: 
	      e->Acceptor->Position + e->offset1 - COFFSET;
	    p2 = (e->evidence)? 
	      e->Donor->Position + e->offset2: 
	      e->Donor->Position + e->offset2 - COFFSET;
	    /* Get the cDNA for this exon */
	    tmpDNA[0] = '\0';
	    rs[0] = '\0';
   
	    for(j=p1; j <= p2; j++)
	      tmpDNA[j-p1] = s[j];
        
	    tmpDNA[j-p1]='\0';
        
	    ReverseSubSequence(0, j-p1-1, tmpDNA, rs);
	    rs[j-p1]='\0';
   
	    /* Concat the current cDNA after the previous */
	    strcat(cDNA,rs);
   
	    *nNN += p2-p1+1;
	 /*  } */
	  e = e->PreviousExon;
		  
	} 
      /* Set free reverse sequence */
      free(rs);
    }

  /* Set free temporary result sequence */
  free(tmpDNA);
}
