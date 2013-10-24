/*************************************************************************
*                                                                        *
*   Module: ComputeTP                                                    *
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


long computeTP (exonGFF* Pred, exonGFF* Real, long nPreds, long nReals)
{
  /* CP is the counter for the predicted exons */
  /* CR is the counter for the real exons */
  /* TP is the number of True Positive nucleotides */
  long CP;
  long CP2;
  long CR;
  long TP;

  /* For each annotation do ... */
  for (  CP = CP2 = CR = TP = 0; CR < nReals && CP < nPreds; CR++)
    {
      /* 1. Skipping pred exons that end before the real exon */ 
      while ((Pred+CP)->Position2 < (Real+CR)->Position1 && CP < nPreds)
	  CP++;

      /* 2. Working with the preds that overlap the real exon */ 
      CP2=CP;
      while ((Pred+CP2)->Position1 <= (Real+CR)->Position2 && CP2 < nPreds)
	{
	  /* Pred exon starts before real exon? */
	  if ((Pred+CP2)->Position1 < (Real+CR)->Position1)
	    {
	      /* Pred exon ends before real exon? */
	      if ((Pred+CP2)->Position2 <= (Real+CR)->Position2)
		TP += (Pred+CP2)->Position2 - (Real+CR)->Position1 + 1;
	      else
		TP += (Real+CR)->Position2 - (Real+CR)->Position1 + 1;
	    }
	  else
	    {
	      if ( (Pred+CP2)->Position2 <= (Real+CR)->Position2 )
		TP += (Pred+CP2)->Position2 - (Pred+CP2)->Position1 + 1;
	      else
		TP += (Real+CR)->Position2 - (Pred+CP2)->Position1 + 1;      
	    }
	  CP2++;   
	}
    }
  return(TP);
}

/** Total stats: nucleotide level **/
void computeTotalNucleotideLevelValues(Svalues* stats)
{
  double RN,PN;
    
  /* temporary vars */
  RN = stats->LengthSequence - stats->CDSreal;
  PN = stats->LengthSequence - stats->CDSpred;

  /* */
  stats->FP = stats->CDSpred - stats->TP;
  stats->FN = stats->CDSreal - stats->TP;
  stats->TN = RN - stats->FP;
  stats->CC=
    (((double)stats->TP * (double)stats->TN) - ((double)stats->FN * (double)stats->FP)) / sqrt((double)stats->CDSreal * (double)RN * (double)stats->CDSpred * (double)PN);
  
  stats->AC=
    (((double)stats->TP / (double)stats->CDSreal  + (double)stats->TP / (double)stats->CDSpred + (double)stats->TN / RN + (double)stats->TN / PN) / 
     4 -0.5) * 2;

  stats->SN= (double)stats->TP / (double)stats->CDSreal;                  /* CDSreal=TP+FN */
  stats->SP= (double) stats->TP / (double)stats->CDSpred;                 /* CDSpred=TP+FP */
} 


/*********   Zentral routine:                                 **************/

 
void computeNucleotideLevelValues(packExons* allExons, Svalues* stats)
{
  double RN,PN;
    
  stats->TP = computeTP(allExons->Pred->exon[FORWARD], 
			allExons->Real->exon[FORWARD], 
			allExons->Pred->numExons[FORWARD], 
			allExons->Real->numExons[FORWARD]);
  
  stats->TP += computeTP(allExons->Pred->exon[REVERSE], 
			allExons->Real->exon[REVERSE], 
			allExons->Pred->numExons[REVERSE], 
			allExons->Real->numExons[REVERSE]);

  /* temporary vars */
  RN = stats->LengthSequence - stats->CDSreal;
  PN = stats->LengthSequence - stats->CDSpred;

  /* */
  stats->FP = stats->CDSpred - stats->TP;
  stats->FN = stats->CDSreal - stats->TP;
  stats->TN = RN - stats->FP;
  stats->CC=
    (((double)stats->TP * (double)stats->TN) - ((double)stats->FN * (double)stats->FP)) / sqrt((double)stats->CDSreal * (double)RN * (double)stats->CDSpred * (double)PN);
  
  stats->AC=
    (((double)stats->TP / (double)stats->CDSreal  + (double)stats->TP / (double)stats->CDSpred + (double)stats->TN / RN + (double)stats->TN / PN) / 
     4 -0.5) * 2;

  stats->SN= (double)stats->TP / (double)stats->CDSreal;                  /* CDSreal=TP+FN */
  stats->SP= (double) stats->TP / (double)stats->CDSpred;                 /* CDSpred=TP+FP */
}



