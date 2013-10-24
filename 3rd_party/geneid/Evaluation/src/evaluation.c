/*************************************************************************
*                                                                        *
*   Module: evaluation                                                   *
*                                                                        *
*   Main program. Management of the actions of evaluation.               *
*                                                                        *
*   This file is part of the evaluation Distribution                     *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Genis   PARRA FARRE                           *
*                          Roderic GUIGO SERRA                           * 
*     with contributions from:                                           *
*                          Moises  BURSET ALVAREDA                       *
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

#include "evaluation.h"

/* Setup Flags */
int   VRB=0;  /* verbose */
int   TST=0;  /* Show total stats (for n>1 sequences) */
int   AST=0;  /* Show average stats (for n>1 sequences) */
int   SST=0;  /* Show a short output only some values */


/*****************************************************************************/

int main (int argc, char *argv[])
{
  /* Input Files */
  char   RealFile[FILENAMELENGTH],
         PredictionsFile[FILENAMELENGTH];
  
  char mess[MAXSTRING];
  
  /* Predictions and annotations */
  packExons* allExons;

  /* Genes information */
  packGenes* allGenes;
  
  /* Statistics values */
  Svalues* stats;
  Svalues* totalStats;

  long nSequences;

  FILE* fPredictionsFile;
  FILE* fRealFile;
  
  int i;

  /* 0. Read setup options */
  readargv(argc,argv,PredictionsFile,RealFile);
  printMess("reading args...");

  printMess("Allocating memory");
  /* 1. Allocating data structures */
  
  if ((allExons = (packExons*) malloc(sizeof(packExons))) == NULL)
    printError("Not enough space to hold exons information"); 

  if ((allExons->Pred = (exons*) malloc(sizeof(exons))) == NULL)
    printError("Not enough space to hold predictions information"); 
  
  if ((allExons->Real = (exons*) malloc(sizeof(exons))) == NULL)
    printError("Not enough space to hold annotations information"); 

  if ((allGenes = (packGenes*) malloc(sizeof(packGenes))) == NULL)
    printError("Not enough space to hold gene information"); 
  
  if ((allGenes->Pred = (genes*) malloc(sizeof(genes))) == NULL)
    printError("Not enough space to hold predicted genes"); 

  if ((allGenes->Real = (genes*) malloc(sizeof(genes))) == NULL)
    printError("Not enough space to hold annotations information"); 
  
  for (i=0; i<STRANDS; i++)
    {
      if ((allExons->Pred->exon[i] = (exonGFF*) calloc(NUMEXONS,sizeof(exonGFF))) == NULL)
	printError("Not enough space to hold predictions information");       

      if ((allExons->Real->exon[i] = (exonGFF*) calloc(NUMEXONS,sizeof(exonGFF))) == NULL)
	printError("Not enough space to hold real information");
      
      if ((allGenes->Pred->gen[i] = (gene*) calloc(MAXGENES,sizeof(gene))) == NULL)
	printError("Not enough space to hold predicted genes");

      if ((allGenes->Real->gen[i] = (gene*) calloc(MAXGENES,sizeof(gene))) == NULL)
	printError("Not enough space to hold real genes");
    }
  
  if ((stats = (Svalues*) malloc(sizeof(Svalues))) == NULL)
    printError("Not enough space to hold Svalues"); 

  if ((totalStats = (Svalues*) malloc(sizeof(Svalues))) == NULL)
    printError("Not enough space to hold total Svalues"); 
  
  
  printMess("**** Running evaluation 1.0.  2001 by Enrique Blanco (www1.imim.es)");
 
  /* Open both files */
  fPredictionsFile = OpenFile(PredictionsFile);
  fRealFile = OpenFile(RealFile);

  for(nSequences = 0; !feof(fPredictionsFile) && !feof(fRealFile); nSequences++)
    {
      /* 2. Extract information about current sequence */
      sprintf(mess,"Get info about sequence %ld",nSequences);
      printMess(mess);
      ExtractInfo (fRealFile, stats);

      strcpy(totalStats->Locus,stats->Locus);
  
      sprintf(mess,"Locus: %s -- Length: %ld",stats->Locus,stats->LengthSequence);
      printMess(mess);
  
      /* 3. Reading predictions */
      stats->ExP = ReadExonsGFF(stats->Locus,
				fPredictionsFile,
				allExons->Pred,
				&(stats->CDSpred),
				allGenes->Pred);
      
      sprintf(mess,"Predictions = %ld [%ld + %ld] (%ld genes: %ld + %ld)",
	      stats->ExP,
	      allExons->Pred->numExons[FORWARD],
	      allExons->Pred->numExons[REVERSE],
	      allGenes->Pred->numGenes[FORWARD] + allGenes->Pred->numGenes[REVERSE],
	      allGenes->Pred->numGenes[FORWARD],allGenes->Pred->numGenes[REVERSE]);
      printMess(mess);  
      
      /* 4. Reading annotations */
      stats->ExR = ReadExonsGFF(stats->Locus,
				fRealFile,
				allExons->Real,
				&(stats->CDSreal),
				allGenes->Real);
      sprintf(mess,"Annotations = %ld [%ld + %ld] (%ld genes: %ld + %ld)",
	      stats->ExR,
	      allExons->Real->numExons[FORWARD],
	      allExons->Real->numExons[REVERSE],
	      allGenes->Real->numGenes[FORWARD] + allGenes->Real->numGenes[REVERSE],
	      allGenes->Real->numGenes[FORWARD],allGenes->Real->numGenes[REVERSE]);
      printMess(mess); 

      /* 5. ComputeTP and other Nucleotide Level values */
      printMess("Computing Nucleotide Level statistics...");
      computeNucleotideLevelValues(allExons,stats); 
      
      /* 6. ComputeTPe and other Exon Level values */
      printMess("Computing Exon Level statistics...");
      computeExonLevelValues(allExons,stats);   
      
      /* 7. Compute TPg and other gene level values */
      printMess("Computing gene Level statistics...");
      computeGeneLevelValues(allExons, allGenes,stats);
      
      /* 8. Updating total results */
      updateTotal(stats, totalStats);

      /* 9. Output results */
      Output(stats); 
    }

  FinalOutput(totalStats,nSequences);

  /* 9. The End */
  exit(0);
  return(0);
}


