/*************************************************************************
*                                                                        *
*   Module: CookingGenes                                                 *
*                                                                        *
*   Processing best gene to print by using the selected format           *
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

/*  $Id: CookingGenes.c,v 1.31 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

extern int X10;
extern int GFF;
extern int GFF3;
extern int XML;
extern int cDNA;
extern int PSEQ;
extern int tDNA;
extern int PRINTINT;

/* Local data structure to record stats about every gene */
typedef struct s_gene
{
  long nexons;
  long nintrons;
  long nutrs;
  long ncds;
  long nfeats;
  float score;
  exonGFF *start;
  exonGFF *end;
} gene;

/* Printing a protein/DNA sequence (fasta format) */
void printProt(char* Name,
               long ngen,
               char* prot,
               long nAA,
               int mode,
	       char* GenePrefix)
{
  long j;
  char header[MAXLINE];
  
  /* 1. Print the header(geneid format): protein or genomic sequence */
  if (GFF3){
    if (mode == PROT){
      sprintf(header,"\n>%s_predicted_protein_%s%s_%ld\n",
	      VERSION,
	      GenePrefix,
	      Name,
	      ngen);
    }else{
      
      if (mode == cDNA){
	sprintf(header,"\n>%s_predicted_CDS_%s%s_%ld\n",
	      VERSION,
	      GenePrefix,
	      Name,
	      ngen);
      }else{
	sprintf(header,"\n>%s_predicted_transcript_%s%s_%ld\n",
	      VERSION,
	      GenePrefix,
	      Name,
	      ngen);
      }
    }

  } else {
      if (mode == PROT){
	sprintf(header,"\n>%s%s_%ld|%s_predicted_protein_%ld|%ld_AA\n",
	      GenePrefix,
	      Name,
	      ngen,
	      VERSION,
	      ngen,
	      nAA);
      }else{
	if (mode == cDNA){
	sprintf(header,"\n>%s%s_%ld|%s_predicted_cDNA_%ld|%ld_NN\n",
	      GenePrefix,
	      Name,
	      ngen,
	      VERSION,
	      ngen,
	      nAA);
	}else{
	  sprintf(header,"\n>%s%s_%ld|%s_predicted_transcript_%ld|%ld_NN\n",
	      GenePrefix,
	      Name,
	      ngen,
	      VERSION,
	      ngen,
	      nAA);
	}
      }
  }
  /* Header left out in XML format */
  if (!XML)
    printf("%s",header);
  else
    printf("\t");
  
  /* 2. Print the input sequence */
  for(j=0; j < strlen(prot); j++)
    {
      printf("%c",prot[j]);
      if (!((j+1) % FASTALINE))
	{
	  printf("\n");
	  if (XML)
            printf("\t");
	}
    }
  printf("\n");
  
  /* One while line between 2 genes */
  if (!GFF3 && !XML && (mode == PROT))
    printf("\n");  	
}

/* Returns signal types and profiles according to the type of input exon */
void selectFeatures(char* exonType,
                    char exonStrand,
                    profile** p1,
                    profile** p2,
                    int* type1,
                    int* type2,
                    int* strand,
                    gparam* gp)
{
  if (exonStrand == '+')
    {
      *strand = FORWARD;
      if (!strcmp(exonType,sFIRST))
	{
	  *type1 = STA;
	  *type2 = DON;
	  *p1 = gp->StartProfile;
	  *p2 = gp->DonorProfile;
	}
      if (!strcmp(exonType,sINTERNAL))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;    
	}
      if (!strcmp(exonType,sZEROLENGTH))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;    
	}
      if (!strcmp(exonType,sTERMINAL))
	{
	  *type1 = ACC;
	  *type2 = STO;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->StopProfile;    
	}
      if (!strcmp(exonType,sSINGLE))
	{
	  *type1 = STA;
	  *type2 = STO;
	  *p1 = gp->StartProfile;
	  *p2 = gp->StopProfile;    
	}
      if (!strcmp(exonType,sINTRON))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRINTRON))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTR5INTRON))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTR3INTRON))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRFIRST))
	{
	  *type1 = TSS;
	  *type2 = DON;
	  *p1 = gp->StartProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRFIRSTHALF))
	{
	  *type1 = TSS;
	  *type2 = STA;
	  *p1 = gp->StartProfile;
	  *p2 = gp->StartProfile;		      
	}
      if (!strcmp(exonType,sUTRINTERNAL))
	{
	  *type1 = ACC;
	  *type2 = DON;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTR5INTERNALHALF))
	{
	  *type1 = ACC;
	  *type2 = STA;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->StartProfile;		      
	}
      if (!strcmp(exonType,sUTR3INTERNALHALF))
	{
	  *type1 = STO;
	  *type2 = DON;
	  *p1 = gp->StopProfile;
	  *p2 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRTERMINAL))
	{
	  *type1 = ACC;
	  *type2 = TES;
	  *p1 = gp->AcceptorProfile;
	  *p2 = gp->StopProfile;		      
	}
      if (!strcmp(exonType,sUTRTERMINALHALF))
	{
	  *type1 = STO;
	  *type2 = TES;
	  *p1 = gp->StopProfile;
	  *p2 = gp->StopProfile;		      
	}
    }
  /* Reverse strand */
  else
    {
      *strand = REVERSE;
      if (!strcmp(exonType,sFIRST))
	{
	  *type2 = STA;
	  *type1 = DON;
	  *p2 = gp->StartProfile;
	  *p1 = gp->DonorProfile;    
	}
      if (!strcmp(exonType,sINTERNAL))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;    
	}
      if (!strcmp(exonType,sZEROLENGTH))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;    
	}
      if (!strcmp(exonType,sINTRON))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;
	}
      if (!strcmp(exonType,sUTRINTRON))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;
	}
      if (!strcmp(exonType,sUTR5INTRON))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;
	}
      if (!strcmp(exonType,sUTR3INTRON))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;
	}
      if (!strcmp(exonType,sTERMINAL))
	{
	  *type2 = ACC;
	  *type1 = STO;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->StopProfile;    
	}
      if (!strcmp(exonType,sSINGLE))
	{
	  *type2 = STA;
	  *type1 = STO;
	  *p2 = gp->StartProfile;
	  *p1 = gp->StopProfile;    
	}
       if (!strcmp(exonType,sUTRFIRST))
	{
	  *type2 = TSS;
	  *type1 = DON;
	  *p2 = gp->StartProfile;
	  *p1 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRFIRSTHALF))
	{
	  *type2 = TSS;
	  *type1 = STA;
	  *p2 = gp->StartProfile;
	  *p1 = gp->StartProfile;		      
	}
      if (!strcmp(exonType,sUTRINTERNAL))
	{
	  *type2 = ACC;
	  *type1 = DON;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTR5INTERNALHALF))
	{
	  *type2 = ACC;
	  *type1 = STA;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->StartProfile;		      
	}
      if (!strcmp(exonType,sUTR3INTERNALHALF))
	{
	  *type2 = STO;
	  *type1 = DON;
	  *p2 = gp->StopProfile;
	  *p1 = gp->DonorProfile;		      
	}
      if (!strcmp(exonType,sUTRTERMINAL))
	{
	  *type2 = ACC;
	  *type1 = TES;
	  *p2 = gp->AcceptorProfile;
	  *p1 = gp->StopProfile;		      
	}
      if (!strcmp(exonType,sUTRTERMINALHALF))
	{
	  *type2 = STO;
	  *type1 = TES;
	  *p2 = gp->StopProfile;
	  *p1 = gp->StopProfile;		      
	}
    }
}

/* Processing of genes to make pretty printing afterwards */
/* artScore is the value that must be substracted from total score (forced evidences) */
long CookingInfo(exonGFF* eorig,
                 gene info[],
                 double* artScore)
{
  /* Identifier of current gene */
  long igen;
  int stop,stop1,stop2,currentIsUTR;
  exonGFF* e;
/*   char mess[MAXSTRING]; */
  /* Reset counters into the gene information structure */
  currentIsUTR=0;
  *artScore = 0.0;
  e = eorig;
  for(igen=0; igen < MAXGENE; igen++)
    {
      info[igen].nexons = 0;
      info[igen].nutrs = 0;
      info[igen].ncds = 0;
      info[igen].nintrons = 0;
      info[igen].nfeats = 0;
      info[igen].score = 0.0;
    }
  
  /* Pointer jumping back travelling across exons of multiple genes */
  /* starting from the last exon of the last gene (bottom-up) */
  igen = 0;
  stop = (e->Strand == '*');
  while (!stop)
    {
      /* A. Skip BEGIN/END features */
      if (!strcmp(e->Type,sEND) || !strcmp(e->Type,sBEGIN))
	{
	  /* Skip this feature: substract the score from the total score */
	  *artScore = *artScore + MAXSCORE;

	  /* JUMP! */
	  e = (e-> PreviousExon);
	  stop = (e->Strand == '*');
	}
      else
	{	    
	  /* C. Reverse Genes: (BOTTOM) First->Internal->...->Terminal (TOP) */
	  if (e->Strand == '-')
	    {	
	      info[igen].start = e;
	      info[igen].end = e;
	      info[igen].nfeats++;
	      if (strcmp(e->Type,sINTRON)
		  && strcmp(e->Type,sUTRINTRON)
		  && strcmp(e->Type,sUTR5INTRON)
		  && strcmp(e->Type,sUTR3INTRON)
		  && strcmp(e->Type,sUTRFIRSTHALF)
		  && strcmp(e->Type,sUTR5INTERNALHALF)
		  && strcmp(e->Type,sUTR3INTERNALHALF)
		  && strcmp(e->Type,sUTRTERMINALHALF)){
		info[igen].nexons++;
	      }
	      if (strcmp(e->Type,sINTRON)
		  && strcmp(e->Type,sUTRINTRON)
		  && strcmp(e->Type,sUTR5INTRON)
		  && strcmp(e->Type,sUTR3INTRON)
		  && strcmp(e->Type,sFIRST)
		  && strcmp(e->Type,sINTERNAL)
		  && strcmp(e->Type,sTERMINAL)
		  && strcmp(e->Type,sSINGLE)){
		info[igen].nutrs++;
	      }

	      if (!strcmp(e->Type,sFIRST)
		  || !strcmp(e->Type,sINTERNAL)
		  || !strcmp(e->Type,sTERMINAL)
		  || !strcmp(e->Type,sSINGLE)
		  ){
		info[igen].ncds++;
	      }
	      /* Evidences (annotations) not added if infinitum score */
	      if (e->Score==MAXSCORE)
		*artScore = *artScore + MAXSCORE;
	      else
		info[igen].score += e->Score;
	
	      /* if(!strcmp(e->Type,sUTRFIRSTHALF)||!strcmp(e->Type,sUTR5INTERNALHALF)){ */
		if(!strcmp(e->Type,sUTRFIRSTHALF)||!strcmp(e->Type,sUTR5INTERNALHALF)||!strcmp(e->Type,sUTRFIRST)){
		currentIsUTR = 1;
				    
	      }else{
		currentIsUTR = 0;
	      } 			  
	      /* JUMP! */
	      e = (e-> PreviousExon);       
	      /* stop means end of processing */
	      stop = (e->Strand == '*');
	      /* stop1 means change of gene: new gene found */
	      stop1 = ((!strcmp(e->Type,sFIRST)&&!currentIsUTR) ||
		       (!strcmp(e->Type,sSINGLE)&&!currentIsUTR) || 
		       (!strcmp(e->Type,sUTR5INTERNALHALF)&&!currentIsUTR) || 
		       !strcmp(e->Type,sPROMOTER) || 
		       !strcmp(e->Type,sBEGIN) ||
		       !strcmp(e->Type,sUTRFIRST) ||
		       !strcmp(e->Type,sUTRFIRSTHALF) ||
		       !strcmp(e->Type,sBEGIN) ||
		       e->Strand == '+');
	      while( !stop && !stop1 )
				    
		{  
		  info[igen].nfeats++;
		  if (strcmp(e->Type,sINTRON)
		      && strcmp(e->Type,sUTRINTRON)
		      && strcmp(e->Type,sUTR5INTRON)
		      && strcmp(e->Type,sUTR3INTRON)
		      && strcmp(e->Type,sUTRFIRSTHALF)
		      && strcmp(e->Type,sUTR5INTERNALHALF)
		      && strcmp(e->Type,sUTR3INTERNALHALF)
		      && strcmp(e->Type,sUTRTERMINALHALF)){
		    info[igen].nexons++;
		  }
		  if (strcmp(e->Type,sINTRON)
		      && strcmp(e->Type,sUTRINTRON)
		      && strcmp(e->Type,sUTR5INTRON)
		      && strcmp(e->Type,sUTR3INTRON)
		      && strcmp(e->Type,sFIRST)
		      && strcmp(e->Type,sINTERNAL)
		      && strcmp(e->Type,sTERMINAL)
		      && strcmp(e->Type,sSINGLE)){
		    info[igen].nutrs++;
		  }
		  if (!strcmp(e->Type,sFIRST)
		      || !strcmp(e->Type,sINTERNAL)
		      || !strcmp(e->Type,sTERMINAL)
		      || !strcmp(e->Type,sSINGLE)
		      ){
		    info[igen].ncds++;
		  }
		  /* Evidences (annotations) not summed if infinite score */
		  if (e->Score==MAXSCORE)
		    *artScore = *artScore + MAXSCORE;
		  else
		    info[igen].score += e->Score;
		  info[igen].end = e;
		  if(!strcmp(e->Type,sUTRFIRSTHALF)||!strcmp(e->Type,sUTR5INTERNALHALF)||!strcmp(e->Type,sUTRFIRST)){
		/*if(!strcmp(e->Type,sUTRFIRSTHALF)||!strcmp(e->Type,sUTR5INTERNALHALF)){*/
		    currentIsUTR = 1;
				    
		  }else{
		    currentIsUTR = 0;
		  }			  
		  /* JUMP loop! */
		  e = (e-> PreviousExon);
		  stop = (e->Strand == '*');
		  stop1 = ((!strcmp(e->Type,sFIRST)&&!currentIsUTR) ||  
			   (!strcmp(e->Type,sSINGLE)&&!currentIsUTR) ||
			   !strcmp(e->Type,sPROMOTER) || 
			   !strcmp(e->Type,sBEGIN) ||
			   !strcmp(e->Type,sUTRFIRST) ||
			   !strcmp(e->Type,sUTRFIRSTHALF) || 
			   e->Strand == '+'); 
					  
		  
		} 
	    }
	  else
	    {
	      /* D. Forward Genes: (BOTTOM) Terminal->Internal->...->First (TOP) */
	      if (e->Strand == '+')
		{
		  info[igen].start = e;
		  info[igen].end = e;
		  info[igen].nfeats++;
		  if (strcmp(e->Type,sINTRON)
		      && strcmp(e->Type,sUTRINTRON)
		      && strcmp(e->Type,sUTR5INTRON)
		      && strcmp(e->Type,sUTR3INTRON)
		      && strcmp(e->Type,sUTRFIRSTHALF)
		      && strcmp(e->Type,sUTR5INTERNALHALF)
		      && strcmp(e->Type,sUTR3INTERNALHALF)
		      && strcmp(e->Type,sUTRTERMINALHALF)){
		    info[igen].nexons++;
		  }
		  if (strcmp(e->Type,sINTRON)
		      && strcmp(e->Type,sUTRINTRON)
		      && strcmp(e->Type,sUTR5INTRON)
		      && strcmp(e->Type,sUTR3INTRON)
		      && strcmp(e->Type,sFIRST)
		      && strcmp(e->Type,sINTERNAL)
		      && strcmp(e->Type,sTERMINAL)
		      && strcmp(e->Type,sSINGLE)){
		    info[igen].nutrs++;
		  }
		  if (!strcmp(e->Type,sFIRST)
		      || !strcmp(e->Type,sINTERNAL)
		      || !strcmp(e->Type,sTERMINAL)
		      || !strcmp(e->Type,sSINGLE)
		      ){
		    info[igen].ncds++;
		  }
		  if (e->Score==MAXSCORE)
		    *artScore = *artScore + MAXSCORE;
		  else
		    info[igen].score += e->Score;
			
		  if(!strcmp(e->Type,sUTRTERMINALHALF)||!strcmp(e->Type,sUTR3INTERNALHALF)){
		    currentIsUTR = 1;
				    
		  }else{
		    currentIsUTR = 0;
		  }
		  /* JUMP */
		  e = (e-> PreviousExon);  
		  stop = (e->Strand == '*');
		  /* stop2 means change of gene */
		  stop2 = ((!strcmp(e->Type,sTERMINAL)&&!currentIsUTR) ||  
			   (!strcmp(e->Type,sSINGLE)&&!currentIsUTR) ||
			   !strcmp(e->Type,sPROMOTER) ||
			   !strcmp(e->Type,sBEGIN) ||
			   !strcmp(e->Type,sUTRTERMINAL) ||
			   !strcmp(e->Type,sUTRTERMINALHALF) ||
			   e->Strand == '-'); 
		
		
		  while( !stop && !stop2 )
		    { 
		      info[igen].nfeats++;
		      if (strcmp(e->Type,sINTRON)
			  && strcmp(e->Type,sUTRINTRON)
			  && strcmp(e->Type,sUTR5INTRON)
			  && strcmp(e->Type,sUTR3INTRON)
			  && strcmp(e->Type,sUTRFIRSTHALF)
			  && strcmp(e->Type,sUTR5INTERNALHALF)
			  && strcmp(e->Type,sUTR3INTERNALHALF)
			  && strcmp(e->Type,sUTRTERMINALHALF)){
			info[igen].nexons++;
		      }
		      if (strcmp(e->Type,sINTRON)
			  && strcmp(e->Type,sUTRINTRON)
			  && strcmp(e->Type,sUTR5INTRON)
			  && strcmp(e->Type,sUTR3INTRON)
			  && strcmp(e->Type,sFIRST)
			  && strcmp(e->Type,sINTERNAL)
			  && strcmp(e->Type,sTERMINAL)
			  && strcmp(e->Type,sSINGLE)){
			info[igen].nutrs++;
		      }
		      if (!strcmp(e->Type,sFIRST)
			  || !strcmp(e->Type,sINTERNAL)
			  || !strcmp(e->Type,sTERMINAL)
			  || !strcmp(e->Type,sSINGLE)
			  ){
			info[igen].ncds++;
		      }
		      /* Evidences (annotations) not added if infinitum score */
		      if (e->Score==MAXSCORE)
			*artScore = *artScore + MAXSCORE;
		      else
			info[igen].score += e->Score;
		      info[igen].end = e;
		      if(!strcmp(e->Type,sUTRTERMINALHALF)||!strcmp(e->Type,sUTR3INTERNALHALF)){
			currentIsUTR = 1;
				    
		      }else{
			currentIsUTR = 0;
		      }				
		      /* JUMP loop! */
		      e = (e-> PreviousExon);
		      stop = (e->Strand == '*');
		      stop2 = ((!strcmp(e->Type,sTERMINAL)&&!currentIsUTR) ||
			       (!strcmp(e->Type,sSINGLE)&&!currentIsUTR) ||
			       !strcmp(e->Type,sPROMOTER) || 
			       !strcmp(e->Type,sBEGIN) ||
			       !strcmp(e->Type,sUTRTERMINAL) ||
			       !strcmp(e->Type,sUTRTERMINALHALF) ||
			       e->Strand == '-'); 
		    }	
		}
	    }
	  info[igen].nintrons =  info[igen].nexons - 1;
	  igen++;
	} 
    }
  
  return (igen);
}

/* Print a gene according to formatted output selected and info structure */
/* taa[exon][0] means first amino acid id. and taa[exon][1] means last one */
void PrintGene(exonGFF* start,
               exonGFF* end,
               exonGFF* last,
               char Name[],
               char* s,
               gparam* gp,
               dict* dAA,
               long igen,
               long nAA,
               int** tAA,
               long cexons,
	       long cintrons,
	       long ccds,
	       long cutrs,
	       long cfeats,
	       gene* info, long geneindex,
	       char* GenePrefix)
{
  exonGFF* eaux;
  profile* p1;
  profile* p2;
  int type1;
  int type2;
  int strand;
  int intronnum = 0;
  int exonnum = 1;
  int cdsnum = 1;
  int utrnum = 1;
  int featnum =0;
/*   char mess[MAXSTRING]; */
/*  if (start != end){ */
   eaux = start->PreviousExon;
/*    if (((!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(last->Type,sUTR3INTERNALHALF)||!strcmp(last->Type,sUTRTERMINALHALF)))){ */
/*      sprintf(mess,"cexons: %ld",cexons);printMess(mess);  */
/*    } */
   if (strcmp(start->Type,sINTRON)
      && strcmp(start->Type,sUTRINTRON)
      && strcmp(start->Type,sUTR5INTRON)
      && strcmp(start->Type,sUTR3INTRON)

      &&!((!strcmp(start->Type,sUTRFIRSTHALF)||!strcmp(start->Type,sUTR5INTERNALHALF))&&(!strcmp(last->Type,sFIRST)||!strcmp(last->Type,sSINGLE)))
      &&!((!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(last->Type,sUTR3INTERNALHALF)||!strcmp(last->Type,sUTRTERMINALHALF)))
       &&!((!strcmp(last->Type,sUTRFIRSTHALF)||!strcmp(last->Type,sUTR5INTERNALHALF))&&(!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sSINGLE)))
       &&!((!strcmp(last->Type,sTERMINAL)||!strcmp(last->Type,sSINGLE))&&(!strcmp(start->Type,sUTR3INTERNALHALF)||!strcmp(start->Type,sUTRTERMINALHALF)))
      
       ){
    cexons++;
		      
  }	      
/*  }else{ */

/*  } */
  
  if (strcmp(start->Type,sINTRON)
      && strcmp(start->Type,sUTRINTRON)
      && strcmp(start->Type,sUTR5INTRON)
      && strcmp(start->Type,sUTR3INTRON)
      && strcmp(start->Type,sFIRST)
      && strcmp(start->Type,sINTERNAL)
      && strcmp(start->Type,sTERMINAL)
      && strcmp(start->Type,sSINGLE)){
    cutrs++;
  }else{
    if (strcmp(start->Type,sINTRON)
	&& strcmp(start->Type,sUTRINTRON)
	&& strcmp(start->Type,sUTR5INTRON)
	&& strcmp(start->Type,sUTR3INTRON)){
      ccds++;
    }
  }
  cfeats++;
/* a. Recursive case */
  if (start != end)
    {
      eaux = start->PreviousExon;
      /*	if (PRINTINT && strcmp(start->Type,sINTRON) && strcmp(start->Type,sUTRINTRON) && strcmp(start->Type,sUTR5INTRON)&& strcmp(start->Type,sUTR3INTRON)&&*/
	if (strcmp(start->Type,sINTRON) && strcmp(start->Type,sUTRINTRON) && strcmp(start->Type,sUTR5INTRON)&& strcmp(start->Type,sUTR3INTRON)&&
	  !(
	    ((start->Strand == '-')&& (!strcmp(start->Type,sUTRFIRSTHALF)||!strcmp(start->Type,sUTR5INTERNALHALF))&&(!strcmp(eaux->Type,sFIRST)||!strcmp(eaux->Type,sSINGLE)))||
	    ((start->Strand == '+')&& (!strcmp(eaux->Type,sUTRFIRSTHALF)||!strcmp(eaux->Type,sUTR5INTERNALHALF))&&(!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sSINGLE)))||
	    ((start->Strand == '-')&& (!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(eaux->Type,sUTR3INTERNALHALF)||!strcmp(eaux->Type,sUTRTERMINALHALF)))||
	    ((start->Strand == '+')&& (!strcmp(eaux->Type,sTERMINAL)||!strcmp(eaux->Type,sSINGLE))&&(!strcmp(start->Type,sUTR3INTERNALHALF)||!strcmp(start->Type,sUTRTERMINALHALF)))
	    )
	  ){
	cintrons++;
      }
      /* a.1. Recursive call to print before the rest of the gene */
	/*cexons = cintrons; +(cfeats - (ccds +cutrs)) + 1;*/
      PrintGene(eaux,end,start,Name,s,gp,dAA,igen,nAA,tAA,cexons,cintrons,ccds,cutrs,cfeats,info,geneindex,GenePrefix);
/*       cexons = cintrons+(cfeats - ccds -cutrs) + 1; */
      if (strcmp(start->Type,sINTRON)&&strcmp(start->Type,sUTRINTRON)&&strcmp(start->Type,sUTR5INTRON)&&strcmp(start->Type,sUTR3INTRON)){
	/* a.2. printing this exon: XML, extend, gff or geneid format */      
	if (XML)
	  {
	    selectFeatures(start->Type,start->Strand,
			   &p1,&p2,&type1,&type2,&strand,gp);
	    if (!strcmp(end->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
		cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
		PrintXMLExon(start,Name,igen,cdsnum,type1,type2,GenePrefix);
	    }else{
		utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
		PrintXMLExon(start,Name,igen,utrnum,type1,type2,GenePrefix);
	    }
	  }
	else 
	  if (X10)
	    {
	      /* Print both sites of exon Start */
	      selectFeatures(start->Type,start->Strand,
			     &p1,&p2,&type1,&type2,&strand,gp);
	      PrintSite(start->Acceptor,type1,Name,strand,s,p1);
	      if (!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
		cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
		PrintGCDS(start,Name,s,dAA,igen,tAA[cfeats - COFFSET][0],tAA[cfeats - COFFSET][1],nAA,cdsnum,GenePrefix);
	      }else{
		utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
		PrintGUTR(start,Name,s,igen,utrnum,GenePrefix);
	      }
	      PrintSite(start->Donor,type2,Name,strand,s,p2);
	    }
	  else {	    
	    /*if (PRINTINT && (strcmp(eaux->Type,sINTRON)&&strcmp(eaux->Type,sUTRINTRON)&&strcmp(eaux->Type,sUTR5INTRON)&&strcmp(eaux->Type,sUTR3INTRON))&&*/

	    if ((strcmp(eaux->Type,sINTRON)&&strcmp(eaux->Type,sUTRINTRON)&&strcmp(eaux->Type,sUTR5INTRON)&&strcmp(eaux->Type,sUTR3INTRON))&&
		!(
		 ((start->Strand == '-')&& (!strcmp(start->Type,sUTRFIRSTHALF)||!strcmp(start->Type,sUTR5INTERNALHALF))&&(!strcmp(eaux->Type,sFIRST)||!strcmp(eaux->Type,sSINGLE)))||
		 ((start->Strand == '+')&& (!strcmp(eaux->Type,sUTRFIRSTHALF)||!strcmp(eaux->Type,sUTR5INTERNALHALF))&&(!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sSINGLE)))||
		 ((start->Strand == '-')&& (!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(eaux->Type,sUTR3INTERNALHALF)||!strcmp(eaux->Type,sUTRTERMINALHALF)))||
		 ((start->Strand == '+')&& (!strcmp(eaux->Type,sTERMINAL)||!strcmp(eaux->Type,sSINGLE))&&(!strcmp(start->Type,sUTR3INTERNALHALF)||!strcmp(start->Type,sUTRTERMINALHALF)))
		 )
		){
	      intronnum = (start->Strand == '-')? cintrons: info[geneindex].nintrons - cintrons + 1;
	      if (PRINTINT){PrintGIntron(eaux,start,Name,igen,intronnum,GenePrefix,0,0.0,start->Type);}
	    }
	    /*this is for CDS and UTR features*/
	    if (!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
	      cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
	      PrintGCDS(start,Name,s,dAA,igen,tAA[cfeats - COFFSET][0],tAA[cfeats- COFFSET][1],nAA,cdsnum,GenePrefix);
	    }else{
	      utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
	      PrintGUTR(start,Name,s,igen,utrnum,GenePrefix);
	    }
	    if(
		 ((start->Strand == '-')&& (!strcmp(start->Type,sUTRFIRSTHALF)||!strcmp(start->Type,sUTR5INTERNALHALF))&&(!strcmp(eaux->Type,sFIRST)||!strcmp(eaux->Type,sSINGLE)))||
		 ((start->Strand == '+')&& (!strcmp(eaux->Type,sUTRFIRSTHALF)||!strcmp(eaux->Type,sUTR5INTERNALHALF))&&(!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sSINGLE)))||
		 ((start->Strand == '-')&& (!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(eaux->Type,sUTR3INTERNALHALF)||!strcmp(eaux->Type,sUTRTERMINALHALF)))||
		 ((start->Strand == '+')&& (!strcmp(eaux->Type,sTERMINAL)||!strcmp(eaux->Type,sSINGLE))&&(!strcmp(start->Type,sUTR3INTERNALHALF)||!strcmp(start->Type,sUTRTERMINALHALF)))
	       ){
	      if (!(((start->Strand == '+')&&(!strcmp(start->Type,sSINGLE))&&(!strcmp(last->Type,sUTR3INTERNALHALF)||!strcmp(last->Type,sUTRTERMINALHALF)))
		    ||
		    ((start->Strand == '-')&&(!strcmp(start->Type,sSINGLE))&&(!strcmp(last->Type,sUTR5INTERNALHALF)||!strcmp(last->Type,sUTRFIRSTHALF))))
		  ){
		if(((start->Strand == '+')&&(!strcmp(eaux->Type,sSINGLE))&&(eaux!=end)&&(!strcmp(eaux->PreviousExon->Type,sUTRFIRSTHALF)||!strcmp(eaux->PreviousExon->Type,sUTR5INTERNALHALF)))
		   ||
		   ((start->Strand == '-')&&(!strcmp(eaux->Type,sSINGLE))&&(eaux!=end)&&(!strcmp(eaux->PreviousExon->Type,sUTRTERMINALHALF)||!strcmp(eaux->PreviousExon->Type,sUTR3INTERNALHALF)))
		   ){
		  
/* 		  sprintf(mess,"cexons: %ld",cexons);printMess(mess); */
		  exonnum = (start->Strand == '-')? (cexons?cexons:1): info[geneindex].nexons - cexons + 1;
/* 		  exonnum = intronnum + 1; */
		  PrintGExon(start,3,Name,igen,exonnum,GenePrefix,0,0.0);
		}else{
		  
/* 		  sprintf(mess,"cexons: %ld",cexons);printMess(mess); */
		  exonnum = (start->Strand == '-')? (cexons?cexons:1): info[geneindex].nexons - cexons + 1;
/* 		  exonnum = intronnum + 1; */
		  PrintGExon(start,2,Name,igen,exonnum,GenePrefix,0,0.0);
		}
	      }
	      
	    }else{
	      if(
		 !(((last->Strand == '-')&& (!strcmp(last->Type,sUTRFIRSTHALF)||!strcmp(last->Type,sUTR5INTERNALHALF))&&(!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sSINGLE)))||
		 ((last->Strand == '+')&& (!strcmp(start->Type,sUTRFIRSTHALF)||!strcmp(start->Type,sUTR5INTERNALHALF))&&(!strcmp(last->Type,sFIRST)||!strcmp(last->Type,sSINGLE)))||
		 ((last->Strand == '-')&& (!strcmp(last->Type,sTERMINAL)||!strcmp(last->Type,sSINGLE))&&(!strcmp(start->Type,sUTR3INTERNALHALF)||!strcmp(start->Type,sUTRTERMINALHALF)))||
		 ((last->Strand == '+')&& (!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE))&&(!strcmp(last->Type,sUTR3INTERNALHALF)||!strcmp(last->Type,sUTRTERMINALHALF)))
		 )){
 		exonnum = (start->Strand == '-')? (cexons?cexons:1): info[geneindex].nexons - cexons + 1; 
		/* exonnum = intronnum + 1; */  
		PrintGExon(start,1,Name,igen,exonnum,GenePrefix,0,0.0);
	      }
	    }
	  }
      }else{
	intronnum = (start->Strand == '-')? cintrons: info[geneindex].nintrons - cintrons + 1;  
	if (PRINTINT){
	  PrintGIntron(eaux,last,Name,igen,intronnum,GenePrefix,1,start->Score,start->Type);	  
	}
      }
    }
  else
    {
      if (strcmp(start->Type,sINTRON)&&strcmp(start->Type,sUTRINTRON)&&strcmp(start->Type,sUTR5INTRON)&&strcmp(start->Type,sUTR3INTRON)){
	/* b. Trivial case: not recursive */
	/* b.1. printing this exon: XML, extend, gff or geneid format */      
	if (XML)
	  {
	    selectFeatures(end->Type,end->Strand,
			   &p1,&p2,&type1,&type2,&strand,gp);
	    if (!strcmp(end->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
		cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
		PrintXMLExon(end,Name,igen,cdsnum,type1,type2,GenePrefix);
	    }else{
		utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
		PrintXMLExon(end,Name,igen,utrnum,type1,type2,GenePrefix);
	    }
	  }
	else 
	  if (X10)
	    {
	      /* Print both sites of exon End */
	      selectFeatures(end->Type,end->Strand,
			     &p1,&p2,&type1,&type2,&strand,gp);
	      PrintSite(end->Acceptor,type1,Name,strand,s,p1);
	      if (!strcmp(end->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
		cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
		featnum = (start->Strand == '-')? cfeats: info[geneindex].nfeats - cfeats + 1;
		PrintGCDS(end,Name,s,dAA,igen,tAA[cfeats - COFFSET][0],tAA[cfeats - COFFSET][1],nAA,cdsnum,GenePrefix);
	      }else{
		utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
		PrintGUTR(end,Name,s,igen,utrnum,GenePrefix);
	      }
	      PrintSite(end->Donor,type2,Name,strand,s,p2);
	    }
	  else {

	    if (!strcmp(end->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){		
	      cdsnum = (start->Strand == '-')? ccds: info[geneindex].ncds - ccds + 1;
	      PrintGCDS(end,Name,s,dAA,igen,tAA[cfeats - COFFSET][0],tAA[cfeats - COFFSET][1],nAA,cdsnum,GenePrefix);
	    }else{
	      utrnum = (start->Strand == '-')? cutrs: info[geneindex].nutrs - cutrs + 1;
	      PrintGUTR(end,Name,s,igen,utrnum,GenePrefix);
	    }
	    if (strcmp(last->Type,sUTRTERMINALHALF)&&strcmp(last->Type,sUTR3INTERNALHALF)&&strcmp(last->Type,sUTRFIRSTHALF)&&strcmp(last->Type,sUTR5INTERNALHALF)
		){
	      if (!strcmp(start->Type,sFIRST)||!strcmp(start->Type,sINTERNAL)||!strcmp(start->Type,sTERMINAL)||!strcmp(start->Type,sSINGLE)){
		exonnum = (start->Strand == '-')? (cexons?cexons:1): info[geneindex].nexons - cexons + 1;
		PrintGExon(end,1,Name,igen,exonnum,GenePrefix,0,0.0);
	      }else{
		if (strcmp(start->Type,sUTRTERMINALHALF)&&strcmp(start->Type,sUTR3INTERNALHALF)&&strcmp(start->Type,sUTRFIRSTHALF)&&strcmp(start->Type,sUTR5INTERNALHALF)
		){
		  exonnum = (start->Strand == '-')? (cexons?cexons:1): info[geneindex].nexons - cexons + 1;
		  PrintGExon(end,1,Name,igen,exonnum,GenePrefix,0,0.0);
		}
	      }
	    }
	  }
      }
    }
} 

/* Main routine: post-processing of predicted genes to display them */
void CookingGenes(exonGFF* e,
                  char Name[],
                  char* s,
                  gparam* gp,
                  dict* dAA,
		  char* GenePrefix)
{
  long igen;
  long ngen;
  gene* info;
  char* prot;
  char* tmpDNA;
  char* tmpTDNA;
  long nAA, nNN, nTN;
  int** tAA;
  double artificialScore;
  long i;
  long cexons;
  long cintrons;
  long cutrs;
  long ccds;
  long cfeats;
/*   char mess[MAXSTRING]; */

  tmpDNA = NULL;
  tmpTDNA = NULL;
  

  /* Get info about each gene */
  if ((info = (gene *) calloc(MAXGENE,sizeof(gene))) == NULL)
    printError("Not enough memory: post-processing genes");
  
  /* tAA[gene][exon[0] is the first amino acid of that exon */
  /* tAA[gene][exon[1] is the last amino acid of that exon */
  /* according to the protein product for that gene */
  if ((tAA = (int**) calloc(MAXEXONGENE,sizeof(int*))) == NULL)
    printError("Not enough memory: tAA general structure");
  
  for(i=0; i<MAXEXONGENE; i++)
    if ((tAA[i] = (int*) calloc(2,sizeof(int))) == NULL)
      printError("Not enough memory: tAA[] structure");
  
  /* Post-processing of genes */
  ngen = CookingInfo(e,info,&artificialScore);
  
  /* Protein space */
  /* if (PSEQ) */
  if ((prot = (char*) calloc(MAXAA,sizeof(char))) == NULL)
    printError("Not enough memory: protein product");
  
  /* cDNA memory if required */
  if (cDNA)
    if ((tmpDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
      printError("Not enough memory: cDNA product");
  /* tDNA memory if required */
  if (tDNA)
    if ((tmpTDNA = (char*) calloc(MAXCDNA,sizeof(char))) == NULL)
      printError("Not enough memory: tDNA product");
  
  /* Principal header: forced annotations not used in gene score sum */
  if (XML)
    printf(" genes=\"%ld\" score =\"%.2f\">\n", 
	   ngen,e -> GeneScore - artificialScore); 
  else
    printf("# Optimal Gene Structure. %ld genes. Score = %.2f \n", 
	   ngen,e -> GeneScore - artificialScore); 
  
  /* Pretty-printing of every gene */
  for(igen=ngen-1; igen>=0; igen--)
    {
      cexons = 0;
      cintrons = 0;
      cutrs = 0;
      ccds = 0;
      cfeats = 0;
      /* Translate gene into protein */
/*       sprintf(mess,"gene: %ld; nfeats: %ld",igen,info[igen].nfeats);printMess(mess); */
      TranslateGene(info[igen].start,s,dAA,(info[igen].nfeats),tAA,prot,&nAA);
      
      /* Get genomic DNA for exons if required */
      if (cDNA)
	GetcDNA(info[igen].start,s,(info[igen].nfeats), tmpDNA, &nNN);
      if (tDNA)
	GetTDNA(info[igen].start,s,(info[igen].nfeats), tmpTDNA, &nTN);
	  
      /* Gene header */
      if (XML)
	printf("   <gene idGene=\"%s%s.G%ld\" strand =\"%s\" exons=\"%ld\" score=\"%.2f\">\n",
	       GenePrefix,
	       Name,
	       ngen-igen,
	       (info[igen].start->Strand == '+')? xmlFORWARD : xmlREVERSE, 
	       info[igen].nexons,
	       info[igen].score);
      else     
	if (strcmp(info[igen].start->Type,sPROMOTER)){
	  printf("# Gene %ld (%s). %ld exons. %ld aa. Score = %.2f \n",
		 ngen-igen,
		 (info[igen].start->Strand == '+')? sFORWARD : sREVERSE,
		 info[igen].nexons,
		 nAA,
		 info[igen].score);
	  if (GFF3){
	    PrintGGene(info[igen].start,info[igen].end,Name,ngen-igen,info[igen].score,GenePrefix);
	    PrintGmRNA(info[igen].start,info[igen].end,Name,ngen-igen,info[igen].score,GenePrefix);
	    if (info[igen].start->Strand == '+'){
	      if (!(!strcmp(info[igen].end->Type,sSINGLE)||!strcmp(info[igen].end->Type,sFIRST)||!strcmp(info[igen].end->Type,sUTRFIRST)||!strcmp(info[igen].end->Type,sUTRFIRSTHALF))){
		/* printf("# 5 prime partial: %s\n",info[igen].end->Type); */
		info[igen].end->five_prime_partial = 1;
	      }
	      if (!(!strcmp(info[igen].start->Type,sSINGLE)||!strcmp(info[igen].start->Type,sTERMINAL)||!strcmp(info[igen].start->Type,sUTRTERMINAL)||!strcmp(info[igen].start->Type,sUTRTERMINALHALF))){
		/* printf("# 3 prime partial: %s\n",info[igen].start->Type); */
		info[igen].start->three_prime_partial = 1;
	      }
	    } else {
	      if (!(!strcmp(info[igen].start->Type,sSINGLE)||!strcmp(info[igen].start->Type,sFIRST)||!strcmp(info[igen].start->Type,sUTRFIRST)||!strcmp(info[igen].start->Type,sUTRFIRSTHALF))){
		/* printf("# 5 prime partial: %s\n",info[igen].start->Type); */
		info[igen].start->five_prime_partial = 1;
	      }
	      if (!(!strcmp(info[igen].end->Type,sSINGLE)||!strcmp(info[igen].end->Type,sTERMINAL)||!strcmp(info[igen].end->Type,sUTRTERMINAL)||!strcmp(info[igen].end->Type,sUTRTERMINALHALF))){
		/* printf("# 3 prime partial: %s\n",info[igen].end->Type); */
		info[igen].end->three_prime_partial = 1;
	      }
	    }
	  } 	 
	} else {
	  printf("# Gene %ld (%s). Promoter. %ld bp\n",
		 ngen-igen,
		 (info[igen].start->Strand == '+')? sFORWARD : sREVERSE,
		 nAA*3);
	}

      PrintGene(info[igen].start, info[igen].end, info[igen].end, Name, s, gp, dAA, ngen-igen,
		nAA,tAA,cexons,cintrons,ccds,cutrs,cfeats,info,igen,GenePrefix);

      if (GFF3)
    	printf ("###\n");
      /* [cDNA] and translated protein */
      if (XML)
	{
	  if (tDNA)
	    {
	      printf("      <tDNA length=\"%ld\">\n",nTN);
	      /* cDNA in FASTA format */
	      printProt(Name,ngen-igen,tmpTDNA,nTN,TDNA,GenePrefix);
	      printf("      </tDNA>\n");
	    }
	  if (cDNA)
	    {
	      printf("      <cDNA length=\"%ld\">\n",nNN);
	      /* cDNA in FASTA format */
	      printProt(Name,ngen-igen,tmpDNA,nNN,DNA,GenePrefix);
	      printf("      </cDNA>\n");
	    }
	  if (PSEQ && nAA > 0) {
	    if (strcmp(info[igen].start->Type,sPROMOTER))
	      {
		printf("      <protein length=\"%ld\">\n",nAA);
		/* Protein in FASTA format */
		printProt(Name,ngen-igen,prot,nAA,PROT,GenePrefix);
		printf("      </protein>\n");
	      }
	  }
	  printf("   </gene>\n");
	}
      else
	if (!(GFF))
	  {
	    /* cDNA */
	    if (tDNA)
	      printProt(Name,ngen-igen,tmpTDNA,nTN,TDNA,GenePrefix);
	    if (cDNA && nNN > 0)
	      printProt(Name,ngen-igen,tmpDNA,nNN,DNA,GenePrefix);
			
	    /* Protein in FASTA format (except promoters) */
	    if (PSEQ && nAA > 0) {
	      if (strcmp(info[igen].start->Type,sPROMOTER))
		printProt(Name,ngen-igen,prot,nAA,PROT,GenePrefix);
	    }
	    else
	      if (!cDNA && !tDNA)
		printf("\n");
	  }
    }

  if (GFF3){
    if (PSEQ || cDNA || tDNA)
      printf("\n##FASTA\n");
    /* Pretty-printing of every gene */
    if (PSEQ) {
      for(igen=ngen-1; igen>=0; igen--)
	{
	  
	  
	  /* Translate gene into protein */
	  TranslateGene(info[igen].start,s,dAA,(info[igen].nfeats),tAA,prot,&nAA);
	  /* Protein in FASTA format (except promoters) */
	  if (strcmp(info[igen].start->Type,sPROMOTER) && nAA>0)
	    printProt(Name,ngen-igen,prot,nAA,PROT,GenePrefix);

	}
    }
    if (cDNA){
      /* Pretty-printing of every gene */
      for(igen=ngen-1; igen>=0; igen--)
	{
	  GetcDNA(info[igen].start,s,(info[igen].nfeats), tmpDNA, &nNN);
	  if (nNN > 0)
	    printProt(Name,ngen-igen,tmpDNA,nNN,DNA,GenePrefix);

	}
    }
    if (tDNA){
      /* Pretty-printing of every gene */
      for(igen=ngen-1; igen>=0; igen--)
	{
	  GetTDNA(info[igen].start,s,(info[igen].nfeats), tmpTDNA, &nTN);
	  printProt(Name,ngen-igen,tmpTDNA,nTN,TDNA,GenePrefix);

	}
    }
  }
	
  /* Freedom operations */
  if (PSEQ)
    free(prot);
  free(info);
  for(i=0; i<MAXEXONGENE; i++)
    free(tAA[i]);
  if (cDNA)
    free(tmpDNA);
}


/*************************************************************************/

/* DEBUG: quick-printing of exons */
void PrintExonGFF (exonGFF *e, char Name[], char Source[])
{
  printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%c\t%hd\t(%s)\n",
         Name,
         Source,
         e->Type,
         e->Acceptor->Position + e->offset1,
         e->Donor->Position + e->offset2,
         e->Score,
         e->Strand,
         e->Frame,
         e->Group);
}

/* DEBUG: quick-printing of genes */
void PrintGeneGFF(exonGFF *e, char Name[], char Source[])
{
  if ((e-> PreviousExon)->Strand != '*') 
    PrintGeneGFF(e->PreviousExon,Name,Source);
  PrintExonGFF (e,Name, Source);
}

