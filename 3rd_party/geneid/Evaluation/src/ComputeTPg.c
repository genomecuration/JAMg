/*************************************************************************
*                                                                        *
*   Module: ComputeTPg                                                   *
*                                                                        *
*   This file is part of the evaluation Distribution                     *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Genis   PARRA  FARRE                          *
*                          Roderic GUIGO  SERRA                          * 
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

long computeTPg (genes* infoPred, genes* infoReal, char* Locus)
{
  long TPg;
  long CP;
  long CR;
  long i;
  int identity;
  int index;

  /* For each strand do ... */
  for (index = TPg = 0; index<STRANDS; index++)
    {
      /* for each predicted gene do ... */
      for (  CP = CR = 0; 
	     CP <  infoPred->numGenes[index] && CR < infoReal->numGenes[index]; CP++)
	{
	  /* Searching the real gene for the prediction */
	  /* Skipping previous gene annotations... */
	  while ((CR < infoReal->numGenes[index]) &&
		 ((infoReal->gen[index] + CR)->Last->Position2 
		  < 
		  (infoPred->gen[index] +CP)->First->Position1))
	    CR++;    
	  
	  /* Verifying if both genes are the same */
	  if ((CR < infoReal->numGenes[index]))
	    if ((infoReal->gen[index]+CR)->numExons == (infoPred->gen[index]+CP)->numExons)
	      {
		identity = TRUE;
		for(i=0; 
		    (i< (infoReal->gen[index]+CR)->numExons) && (identity == TRUE); 
		    i++)
		  {
		    /* Compare each couple of exons */
		    identity  =
		      ((((infoReal->gen[index]+CR)->First)+i)->Position1 == 
		       (((infoPred->gen[index]+CP)->First)+i)->Position1)
		      &&
		      ((((infoReal->gen[index]+CR)->First)+i)->Position2 == 
		       (((infoPred->gen[index]+CP)->First)+i)->Position2);
		  }
		
		if (identity == TRUE)
		  TPg++;
	      }
	}
    }
  return(TPg);
}

long computeMWG(exonGFF* Pred, long nPred, genes* infoReal, int strand)
{
  long MG;
  long CP,CP2;
  long CR;
  long i;
  int missing;
  exonGFF* exon;
  long a,b,lIntersection;

  /* for each real gene do ... */
  for (  CP = CR = MG = 0; 
	 CR <  infoReal->numGenes[strand]; CR++)
    {
      /* for each one of the real exons... */
      for (i=0, exon = (infoReal->gen[strand]+CR)->First, missing=TRUE; 
	   (i<(infoReal->gen[strand] + CR)->numExons) && (missing==TRUE); 
	   i++,exon++)
	{
	  /* Skipping previous predicted exons... */
	  while ((CP < nPred) && (Pred+CP)->Position2 < exon->Position1)
	    CP++;

	  CP2 = CP;
	  /* Select all the possible candidates among the following exons */
	  while ((CP2 < nPred) && (Pred+CP2)->Position1 <= exon->Position2 && (missing==TRUE))
	    {
	      a = MIN(exon->Position2,(Pred+CP2)->Position2);
	      b = MAX(exon->Position1,(Pred+CP2)->Position1);
	      lIntersection = a - b + 1; 

	      /* If they overlap, this is not a Missing Gene */
	      missing = (lIntersection == 0);
	      CP2++;

	    }	 
	}

      if (missing == TRUE)
	MG++;
    }
  return(MG);
}

long computeJSG(genes* infoPred, genes* infoReal)
{
  long frequency;
  long CP,CP2;
  long CR;
  long i,j;
  int overlap;
  exonGFF* exonPred;
  exonGFF* exonReal;
  long a,b,lIntersection;
  int index;

  /* For each strand ... */
  for (index = frequency = 0; index<STRANDS; index++)
    {
      /* for each real gene do ... */
      for (  CP = CR = 0; 
	     CR <  infoReal->numGenes[index] && CP < infoPred->numGenes[index] ; CR++)
	{
	  /* Searching all the genes with one or more exons overlapping with the real gene */
	  /* 1. Skipping previous non-overlapping gene predictions ... */
	  while ((CP < infoPred->numGenes[index]) &&
		 ((infoPred->gen[index]+CP)->Last->Position2 
		  < 
		  (infoReal->gen[index]+CR)->First->Position1))
	    CP++;    
	  
	  CP2 = CP;
	  /* 2. Scan all the possible candidates */
	  while ((CP2 < infoPred->numGenes[index]) &&
		 ((infoPred->gen[index]+CP2)->First->Position1 
		  <= 
		  (infoReal->gen[index]+CR)->Last->Position2))
	    {
	      /* for each one of the predicted exons... */
	      for (i=0, j=0, 
		     exonPred = (infoPred->gen[index]+CP2)->First,
		     exonReal = (infoReal->gen[index]+CR)->First,
		     overlap=FALSE; 
		   (i < (infoPred->gen[index]+CP2)->numExons) && 
		     (j < (infoReal->gen[index]+CR)->numExons) && 
		     (overlap==FALSE); 
		   i++,exonPred++)
		{
		  /* compare predicted exon with every real exon from this gene... */
		  while ((overlap==FALSE) && 
			 (j<(infoReal->gen[index]+CR)->numExons) 
			 && (exonPred->Position2 >= exonReal->Position1))
		    {
		      a = MIN(exonPred->Position2,exonReal->Position2);
		      b = MAX(exonPred->Position1,exonReal->Position1);
		      lIntersection = a - b + 1; 
		      
		      overlap = (lIntersection > 0);
		      
		      /* next real exon */
		      j++;
		      exonReal++;
		    }
		}
	      
	      /* Updating p counter */
	      if (overlap)
		frequency++;
	      
	      /* next predicted gene */
	      CP2++;
	    }
	}
    }

  return(frequency);
}


/** Total stats: exon level **/
void computeTotalGeneLevelValues(Svalues* stats)
{    
  double d1, d2;
                                      
  stats->ratioMG = (double)stats->MG / (double)stats->GeR;

  stats->ratioWG = (double)stats->WG / (double)stats->GeP;

  /* Sensitivity & Sensibility */
  stats->SNg = (double)stats->TPg / (double)stats->GeR;
  stats->SPg = (double)stats->TPg / (double)stats->GeP;
  stats->SNSPg = (stats->SNg + stats->SPg)/2;

  d1 = (double)(stats->GeR - stats->MG);
  d2 = (double)(stats->GeP - stats->WG);
  
  stats->ratioJG = (d2 == 0)? JG_DEFAULT : (double)stats->JG / d2;
  stats->ratioSG = (d1 == 0)? SG_DEFAULT : (double)stats->SG / d1;
}



/*********   Zentral routine:                                 **************/

void computeGeneLevelValues(packExons* allExons,
			    packGenes* allGenes,
			    Svalues* stats)
{    
  double d1, d2;
                                      
  stats->GeR = allGenes->Real->numGenes[FORWARD] + allGenes->Real->numGenes[REVERSE];
  stats->GeP = allGenes->Pred->numGenes[FORWARD] + allGenes->Pred->numGenes[REVERSE];
  
  stats->TPg = computeTPg(allGenes->Pred,allGenes->Real,stats->Locus);

  stats->MG  = computeMWG(allExons->Pred->exon[FORWARD],
			  allExons->Pred->numExons[FORWARD],
			  allGenes->Real,
			  FORWARD);
  stats->MG  += computeMWG(allExons->Pred->exon[REVERSE],
			  allExons->Pred->numExons[REVERSE],
			  allGenes->Real,
			  REVERSE);

  stats->ratioMG = (double)stats->MG / (double)stats->GeR;

  stats->WG  = computeMWG(allExons->Real->exon[FORWARD],
			  allExons->Real->numExons[FORWARD],
			  allGenes->Pred,
			  FORWARD);
  stats->WG  += computeMWG(allExons->Real->exon[REVERSE],
			  allExons->Real->numExons[REVERSE],
			  allGenes->Pred,
			  REVERSE);

  stats->ratioWG = (double)stats->WG / (double)stats->GeP;

  /* Sensitivity & Sensibility */
  stats->SNg = (double)stats->TPg / (double)stats->GeR;
  stats->SPg = (double)stats->TPg / (double)stats->GeP;
  stats->SNSPg = (stats->SNg + stats->SPg)/2;

  d1 = (double)(stats->GeR - stats->MG);
  d2 = (double)(stats->GeP - stats->WG);
 
  stats->JG = computeJSG(allGenes->Real, allGenes->Pred);
  stats->SG = computeJSG(allGenes->Pred, allGenes->Real);
  
  stats->ratioJG = (d2 == 0)? JG_DEFAULT : (double)stats->JG / d2;
  stats->ratioSG = (d1 == 0)? SG_DEFAULT : (double)stats->SG / d1;
}
