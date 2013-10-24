/*************************************************************************
*                                                                        *
*   Module: Output                                                       *
*                                                                        *
*   Management of options for printing results.                          *
*                                                                        *
*   This file is part of the evaluation Distribution                     *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
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

extern int VRB;
extern int TST;
extern int AST;
extern int SST;


/* Printing messages (information) */
void printMess(char* s)
{
  if (VRB)
    fprintf (stderr, "> %s\n",s);
}

/* Printing partial and final results (information) */
void printRes(char* s)
{
  if (VRB)
    fprintf (stderr, "\t%s\n",s);
}

/* Printing error messages */
void printError(char *s)
{
  fprintf(stderr,"Error: %s\n",s);
  exit(1);
}

void updateTotal(Svalues* stats, Svalues* tStats)
{
  /* Number of analyzed bps */
  tStats->LengthSequence += stats->LengthSequence;
  
  /* Nucleotide level */
  tStats->TP += stats->TP;
  tStats->CDSreal += stats->CDSreal;
  tStats->CDSpred += stats->CDSpred;
  
  /* Exon level */
  tStats->ExR += stats->ExR;
  tStats->ExP += stats->ExP;
  tStats->TPe += stats->TPe;
  tStats->ME += stats->ME;
  tStats->WE += stats->WE;
  
  /* Gene level */
  tStats->GeR += stats->GeR;
  tStats->GeP += stats->GeP;
  tStats->TPg += stats->TPg;
  tStats->MG += stats->MG;
  tStats->WG += stats->WG;
  tStats->JG += stats->JG;
  tStats->SG += stats->SG;

  /* Info for computing average */
  tStats->FP += stats->FP;
  tStats->FN += stats->FN;
  tStats->TN += stats->TN;
  tStats->CC += stats->CC;
  tStats->AC += stats->AC;
  tStats->SN += stats->SN;
  tStats->SP += stats->SP;
  
  tStats->SNe += stats->SNe;
  tStats->SPe += stats->SPe;
  tStats->SNSP += stats->SNSP;
  tStats->ratioME += stats->ratioME;
  tStats->ratioWE += stats->ratioWE;
  
  tStats->ratioMG += stats->ratioMG;
  tStats->ratioWG += stats->ratioWG;

  tStats->SNg += stats->SNg;
  tStats->SPg += stats->SPg;
  tStats->SNSPg += stats->SNSPg;
  tStats->ratioJG += stats->ratioJG;
  tStats->ratioSG += stats->ratioSG;
}


/* Printing results ****************/

void OutputNormal(Svalues* stats)
{
  /* Header */
  printf("** Stats for sequence: %8s (%ld bps)\n",
	  stats->Locus,
	  stats->LengthSequence);  
  printf("%s\n",sRULE);

  /* Nucleotide level */
  printf("Nucleotide Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\n",
	 sCDSreal,
	 sCDSpred,
	 sTP,
	 sTN,
	 sFP,
	 sFN,
	 sSN,
	 sSP,
	 sAC,
	 sCC);

  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->CDSreal,
	 stats->CDSpred,
	 stats->TP,
	 stats->TN,	  
	 stats->FP,
	 stats->FN,
	 stats->SN,
	 stats->SP,
	 stats->AC,
	 stats->CC);

  /* Exon level */
  printf("Exon Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sExR,
	 sExP,
	 sTPe,
	 sME,
	 sWE,
	 sratioME,
	 sratioWE,
	 sSNe,
	 sSPe,
	 sSNSP);
  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->ExR,
	 stats->ExP,
	 stats->TPe,
	 stats->ME,
	 stats->WE,
	 stats->ratioME,
	 stats->ratioWE,
	 stats->SNe,
	 stats->SPe,
	 stats->SNSP); 

  /* Gene level */
  printf("Gene Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sGeR,
	 sGeP,
	 sTPg,
	 sMG,
	 sWG,
	 sratioMG,
	 sratioWG,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioJG,
	 sratioSG);

  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->GeR,
	 stats->GeP,
	 stats->TPg,
	 stats->MG,
	 stats->WG,
	 stats->ratioMG,
	 stats->ratioWG,
	 stats->SNg,
	 stats->SPg,
	 stats->SNSPg,
	 stats->ratioJG,
	 stats->ratioSG);
}

void OutputShort(Svalues* stats)
{
  /* Nucleotide level */
  printf("#%10s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sSEQ,
	 sSN,
	 sSP,
	 sCC,
	 sSNe,
	 sSPe,
	 sSNSP,
	 sratioME,
	 sratioWE,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioMG,
	 sratioWG,
	 sratioJG,
	 sratioSG);

  printf(" %10s\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->Locus,
	 stats->SN,
	 stats->SP,
	 stats->CC,
	 stats->SNe,
	 stats->SPe,
	 stats->SNSP,
	 stats->ratioME,
	 stats->ratioWE,
	 stats->SNg,
	 stats->SPg,
	 stats->SNSPg,
	 stats->ratioMG,
	 stats->ratioWG,
	 stats->ratioJG,
	 stats->ratioSG);
}

void OutputShortAverage(Svalues* stats,long nSequences)
{
  printf("#Average:\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sSN,
	 sSP,
	 sCC,
	 sSNe,
	 sSPe,
	 sSNSP,
	 sratioME,
	 sratioWE,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioMG,
	 sratioWG,
	 sratioJG,
	 sratioSG);
  
  printf("\t\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->SN / (double)nSequences,
	 stats->SP / (double)nSequences,
	 stats->CC / (double)nSequences,
	 stats->SNe / (double)nSequences,
	 stats->SPe / (double)nSequences,
	 stats->SNSP / (double)nSequences,
	 stats->ratioME / (double)nSequences,
	 stats->ratioWE / (double)nSequences,
	 stats->SNg / (double)nSequences,
	 stats->SPg / (double)nSequences,
	 stats->SNSPg / (double)nSequences,	 
	 stats->ratioMG / (double)nSequences,
	 stats->ratioWG / (double)nSequences,
	 stats->ratioJG / (double)nSequences,
	 stats->ratioSG / (double)nSequences); 
}

void OutputAverage(Svalues* stats,long nSequences)
{
  /* Header */
  printf("** Average stats for %ld sequences: (%ld analyzed bps)\n",
	 nSequences,
	 stats->LengthSequence);  
  printf("%s\n",sRULE);
  
  /* Nucleotide level */
  printf("Nucleotide Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\n",
	 sCDSreal,
	 sCDSpred,
	 sTP,
	 sTN,
	 sFP,
	 sFN,
	 sSN,
	 sSP,
	 sAC,
	 sCC);
  
  printf(" %6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 (double)stats->CDSreal / (double)nSequences,
	 (double)stats->CDSpred / (double)nSequences,
	 (double)stats->TP / (double)nSequences,
	 (double)stats->TN / (double)nSequences,	  
	 (double)stats->FP / (double)nSequences,
	 (double)stats->FN / (double)nSequences,
	 stats->SN / (double)nSequences,
	 stats->SP / (double)nSequences,
	 stats->AC / (double)nSequences,
	 stats->CC / (double)nSequences);
  
  /* Exon level */
  printf("Exon Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sExR,
	 sExP,
	 sTPe,
	 sME,
	 sWE,
	 sratioME,
	 sratioWE,
	 sSNe,
	 sSPe,
	 sSNSP);
  printf(" %6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 (double)stats->ExR / (double)nSequences,
	 (double)stats->ExP / (double)nSequences,
	 (double)stats->TPe / (double)nSequences,
	 (double)stats->ME / (double)nSequences,
	 (double)stats->WE / (double)nSequences,
	 stats->ratioME / (double)nSequences,
	 stats->ratioWE / (double)nSequences,
	 stats->SNe / (double)nSequences,
	 stats->SPe / (double)nSequences,
	 stats->SNSP / (double)nSequences); 
  
  /* Gene level */
  printf("Gene Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sGeR,
	 sGeP,
	 sTPg,
	 sMG,
	 sWG,
	 sratioMG,
	 sratioWG,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioJG,
	 sratioSG);
  
  printf(" %6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 (double)stats->GeR / (double)nSequences,
	 (double)stats->GeP / (double)nSequences,
	 (double)stats->TPg / (double)nSequences,
	 (double)stats->MG / (double)nSequences,
	 (double)stats->WG / (double)nSequences,
	 stats->ratioMG / (double)nSequences,
	 stats->ratioWG / (double)nSequences,
	 stats->SNg / (double)nSequences,
	 stats->SPg / (double)nSequences,
	 stats->SNSPg / (double)nSequences,
	 stats->ratioJG / (double)nSequences,
	 stats->ratioSG / (double)nSequences); 
}

void OutputShortTotal(Svalues* stats,long nSequences)
{
  printf("#Total:  \t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sSN,
	 sSP,
	 sCC,
	 sSNe,
	 sSPe,
	 sSNSP,
	 sratioME,
	 sratioWE,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioMG,
	 sratioWG,
	 sratioJG,
	 sratioSG);
  
  printf("\t\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->SN,
	 stats->SP,
	 stats->CC,
  	 stats->SNe,
	 stats->SPe,
	 stats->SNSP,
	 stats->ratioME,
	 stats->ratioWE,
	 stats->SNg,
	 stats->SPg,
	 stats->SNSPg,
	 stats->ratioMG,
	 stats->ratioWG,
	 stats->ratioJG,
	 stats->ratioSG); 
}


void OutputTotal(Svalues* stats,long nSequences)
{
  /* Header */
  printf("** TOTAL stats for %ld sequences: (%ld bps)\n",
	 nSequences,
	 stats->LengthSequence);  
  printf("%s\n",sRULE);
  
  /* Nucleotide level */
  printf("Nucleotide Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\n",
	 sCDSreal,
	 sCDSpred,
	 sTP,
	 sTN,
	 sFP,
	 sFN,
	 sSN,
	 sSP,
	 sAC,
	 sCC);
  
  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->CDSreal,
	 stats->CDSpred,
	 stats->TP,
	 stats->TN,	  
	 stats->FP,
	 stats->FN,
	 stats->SN,
	 stats->SP,
	 stats->AC,
	 stats->CC);
  
  /* Exon level */
  printf("Exon Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sExR,
	 sExP,
	 sTPe,
	 sME,
	 sWE,
	 sratioME,
	 sratioWE,
	 sSNe,
	 sSPe,
	 sSNSP);
  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n\n",
	 stats->ExR,
	 stats->ExP,
	 stats->TPe,
	 stats->ME,
	 stats->WE,
	 stats->ratioME,
	 stats->ratioWE,
	 stats->SNe,
	 stats->SPe,
	 stats->SNSP); 
  
  /* Gene level */
  printf("Gene Level\n");
  printf("%s\n",sRULE2);
  printf(" %6s\t%6s\t%6s\t%6s\t%6s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\t%4s\n",
	 sGeR,
	 sGeP,
	 sTPg,
	 sMG,
	 sWG,
	 sratioMG,
	 sratioWG,
	 sSNg,
	 sSPg,
	 sSNSPg,
	 sratioJG,
	 sratioSG);
      
  printf(" %6ld\t%6ld\t%6ld\t%6ld\t%6ld\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n",
	 stats->GeR,
	 stats->GeP,
	 stats->TPg,
	 stats->MG,
	 stats->WG,
	 stats->ratioMG,
	 stats->ratioWG,
	 stats->SNg,
	 stats->SPg,
	 stats->SNSPg,
	 stats->ratioJG,
	 stats->ratioSG); 
}

void Output(Svalues* stats)
{
  if (SST)
    OutputShort(stats);
  else
    OutputNormal(stats);  
}

void FinalOutput(Svalues* stats, long nSequences)
{
  if (AST)
    {
      if (SST)
	OutputShortAverage(stats,nSequences);
      else
	OutputAverage(stats,nSequences);
    }
  if (TST)
    {
      computeTotalNucleotideLevelValues(stats);
      computeTotalExonLevelValues(stats);
      computeTotalGeneLevelValues(stats);
      
      if (SST)
	OutputShortTotal(stats,nSequences);
      else
	OutputTotal(stats,nSequences);
    }
}
