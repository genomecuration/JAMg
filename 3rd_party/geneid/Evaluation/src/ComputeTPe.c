/*************************************************************************
*                                                                        *
*   Module: ComputeTPe                                                   *
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

long computeTPe (exonGFF* Pred, exonGFF* Real, long nPreds, long nReals)
{
  /* CP is the counter for the predicted exons */
  /* CR is the counter for the real exons */
  /* TPe is the number of True Positive exons */
  long CP;
  long CR;
  long TPe;
 
  /* For each prediction do ... */
  for (  CP = CR = TPe = 0; CP <  nPreds && CR < nReals; CP++)
    {
      /* 1. Skipping real exons that end before the current prediction */ 
      while ((Real+CR)->Position2 < (Pred+CP)->Position1 && CR < nReals)
	CR++;

      /* 2. Use the current annotation to measure TPe */
      if (CR < nReals)
	if ( ((Real+CR)->Position1 == (Pred+CP)->Position1)
	     && ((Real+CR)->Position2 == (Pred+CP)->Position2) )
	  TPe++;
      
      /* Reaching end of predictions means STOP */
    }
  
  return(TPe);
}

long computeRatioList(exonGFF* Real, exonGFF* Pred, long nReals, long nPreds)
{
  long CP;
  long CR;
  long nNotFound;
  long a,b,lIntersection;
  
  /* For each annotation do ... */
  for (  CP = CR = nNotFound = 0; CR <  nReals ; CR++)
    {
      /* 1. Skipping pred exons that end before the real exon */ 
      while ((Pred+CP)->Position2 < (Real+CR)->Position1 && CP < nPreds)
          CP++;
  
      /* 2. Testing whether they overlap or not */
      a = MIN((Real+CR)->Position2,(Pred+CP)->Position2);
      b = MAX((Real+CR)->Position1,(Pred+CP)->Position1);
      lIntersection = a - b + 1;
      
      if (lIntersection > 0)
	;
      else
	nNotFound++;
    }
  
  return(nNotFound);
}


/** Total stats: nucleotide level **/
void computeTotalExonLevelValues(Svalues* stats)
{
  stats->ratioME = (double)stats->ME / (double)stats->ExR;

  stats->ratioWE = (double)stats->WE / (double)stats->ExP;

  /* Sensitivity */
  stats->SNe= (double)stats->TPe / (double)stats->ExR;                  

  /* Sensibility */
  stats->SPe= (double) stats->TPe / (double)stats->ExP;                 

  stats->SNSP = (stats->SNe + stats->SPe)/2;
}


/*********   Zentral routine:                                 **************/

void computeExonLevelValues(packExons* allExons, Svalues* stats)
{
  stats->TPe = computeTPe(allExons->Pred->exon[FORWARD],
			  allExons->Real->exon[FORWARD],
			  allExons->Pred->numExons[FORWARD],
			  allExons->Real->numExons[FORWARD]);

  stats->TPe += computeTPe(allExons->Pred->exon[REVERSE],
			  allExons->Real->exon[REVERSE],
			  allExons->Pred->numExons[REVERSE],
			  allExons->Real->numExons[REVERSE]);

  stats->ME = computeRatioList(allExons->Real->exon[FORWARD],
			       allExons->Pred->exon[FORWARD],
			       allExons->Real->numExons[FORWARD],
			       allExons->Pred->numExons[FORWARD]);

  

  stats->ME += computeRatioList(allExons->Real->exon[REVERSE],
				allExons->Pred->exon[REVERSE],
				allExons->Real->numExons[REVERSE],
				allExons->Pred->numExons[REVERSE]);

  stats->ratioME = (double)stats->ME / (double)stats->ExR;

  stats->WE = computeRatioList(allExons->Pred->exon[FORWARD],
			       allExons->Real->exon[FORWARD],
			       allExons->Pred->numExons[FORWARD],
			       allExons->Real->numExons[FORWARD]);

  stats->WE += computeRatioList(allExons->Pred->exon[REVERSE],
				allExons->Real->exon[REVERSE],
				allExons->Pred->numExons[REVERSE],
				allExons->Real->numExons[REVERSE]);

  stats->ratioWE = (double)stats->WE / (double)stats->ExP;

  /* Sensitivity */
  stats->SNe= (double)stats->TPe / (double)stats->ExR;                  

  /* Sensibility */
  stats->SPe= (double) stats->TPe / (double)stats->ExP;                 

  stats->SNSP = (stats->SNe + stats->SPe)/2;
}
