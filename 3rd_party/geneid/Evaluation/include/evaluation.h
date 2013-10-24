/*************************************************************************
*                                                                        *
*   Module: evaluation.h                                                 *
*                                                                        *
*   Main program. Management of the actions of evaluation.               *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Enrique BLANCO GARCIA                         *
*                          Genis   PARRA  FARRE                          *
*                          Roderic GUIGO  SERRA                          *
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

/* Include libraries */
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <malloc.h>  
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <sys/stat.h>

#define LOCUSLENGTH 50                 /* maximum number of chars (locus) */ 
#define FILENAMELENGTH 500             /* maximum length of filenames     */
#define MAXSTRING 500
#define MAXLINE 1000                   /* Max number of chars/inputline   */
#define INFI 9999999                   /* the biggest number of the world */
   
#define NOSCORE -1
#define NOFRAME -1
#define NUMEXONS 1000000
#define NOINFO   0
#define INFO     1
#define NOGROUP  "No_group"
#define MAXGENES 2000
#define NOGENES   -1
#define JG_DEFAULT 1
#define SG_DEFAULT 1

#define STRANDS  2
#define FORWARD  0                     /* DNA - Strands                   */
#define REVERSE  1
  
#define TRUE     1
#define FALSE    0

#define MIN(a,b) (a<b)?a:b;
#define MAX(a,b) (a>b)?a:b;

#define sSEQ      "LocusName"
#define sTP       "TP"
#define sCDSreal  "NuR"
#define sCDSpred  "NuP"
#define sFP       "FP"
#define sFN       "FN"
#define sTN       "TN"
#define sCC       "CC"
#define sAC       "AC"
#define sSN       "SN" 
#define sSP       "SP" 
#define sExR      "ExR"
#define sExP      "ExP"
#define sTPe      "TPe"
#define sSNe      "SNe" 
#define sSPe      "SPe" 
#define sSNSP     "SNSP" 
#define sME       "ME"
#define sWE       "WE"
#define sratioME  "raME"
#define sratioWE  "raWE"
#define sGeR      "GeR"
#define sGeP      "GeP"
#define sTPg      "TPg"
#define sMG       "MG"
#define sWG       "WG"
#define sratioMG  "raMG"
#define sratioWG  "raWG"
#define sSNSPg    "SNSPg"
#define sSNg      "SNg"
#define sSPg      "SPg"
#define sJG       "JG"
#define sSG       "SG"
#define sratioJG  "raJG"
#define sratioSG  "raSG"

#define sRULE     "_________________________________________________\n"
#define sRULE2    "-----------------\n"

typedef struct s_exonGFF
{
  long    Position1;
  long    Position2;
  char    Type[MAXSTRING];
  short   Frame;
  char    Strand;
  double  Score;
} exonGFF;

typedef struct s_exons
{
  exonGFF* exon[STRANDS];
  long numExons[STRANDS];
} exons;

typedef struct s_packExons
{
  exons* Pred;
  exons* Real;
} packExons;


typedef struct s_gene
{
  long numExons;
  exonGFF* First;
  exonGFF* Last;
} gene;

typedef struct s_genes
{
  gene* gen[STRANDS];
  long numGenes[STRANDS];
} genes;

typedef struct s_packGenes
{
  genes* Pred;
  genes* Real;
} packGenes;

typedef struct s_Svalues
{
  char Locus[LOCUSLENGTH];
  long LengthSequence;

  /* Nucleotide level */
  long TP;
  long CDSreal;
  long CDSpred;
  long FP;
  long FN;
  long TN;
  double CC;
  double AC;
  double SN;
  double SP;

  /* Exon level */
  long ExR;
  long ExP;
  long TPe;
  double SNe;
  double SPe;
  double SNSP;
  long ME;
  long WE;
  double ratioME;
  double ratioWE;

  /* Gene level */
  long GeR;
  long GeP;
  long TPg;
  long MG;
  long WG;
  double ratioMG;
  double ratioWG;
  double SNg;
  double SPg;
  double SNSPg;
  long JG;
  long SG;
  double ratioJG;
  double ratioSG;
} Svalues;



/* Headers for evaluation functions */

void readargv (int argc,char* argv[],
               char* PredictionsFile,
               char* RealFile);

long ReadExonsGFF (char* Locus,
		   FILE* file,
                   exons* exons,
                   long* nucleotides,
                   genes* genes);   

void ExtractInfo (FILE* file, Svalues* stats);

FILE* OpenFile(char *FileName);

void printMess(char* s);

void printError(char* s);

void Output(Svalues* stats);   

void FinalOutput(Svalues* stats, long nSequences);   

void computeNucleotideLevelValues(packExons* allExons, Svalues* stats);

void computeExonLevelValues(packExons* allExons, Svalues* stats);

void computeGeneLevelValues(packExons* allExons,
                            packGenes* allGenes,
                            Svalues* stats);

void updateTotal(Svalues* stats, Svalues* tStats);

void computeTotalNucleotideLevelValues(Svalues* stats);

void computeTotalExonLevelValues(Svalues* stats);

void computeTotalGeneLevelValues(Svalues* stats);

