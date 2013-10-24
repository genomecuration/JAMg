/*************************************************************************
*                                                                        *
*   Module: ReadExonsGFF                                                 *
*                                                                        *
*   It reads exons in GFF format.                                        *
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

FILE* OpenFile(char* FileName) 
{
  FILE *file;
  char mess[MAXSTRING];
  
  sprintf(mess,"File %s cannot be open for read",FileName);
  
  if ((file=fopen(FileName, "r"))==NULL)
    printError(mess);

  return(file);
}


void ExtractInfo (FILE *file, Svalues* stats)
{

  char line[MAXLINE];

  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;

  
  line6 = fgets(line,MAXLINE,file);
  while(line6 != NULL && (line[0]=='#' || line[0]=='\n'))
    {
      /* Skip this line */
      printMess("Skipping comment line");
      line6 = fgets(line,MAXLINE,file);
    }

  if (line6!=NULL)
    {
      line1 = (char *) strtok(line,"\t");
      line2 = (char *) strtok(NULL,"\t");
      line3 = (char *) strtok(NULL,"\t"); 
      line4 = (char *) strtok(NULL,"\t");
      line5 = (char *) strtok(NULL,"\t");
      line6 = (char *) strtok(NULL,"\n");
            
      if (sscanf(line1,"%s",stats->Locus) != 1)
		printError("Error on information-about-sequence line: Locusname");
      
      if (sscanf(line5,"%ld",&(stats->LengthSequence)) != 1)
		printError("Error on information-about-sequence line: LengthSequence");
    }
  else
    printError("Error on information-about-sequence line");

}

/* Returns the total number of exons read and the total number of nucleotides */
long ReadExonsGFF (char* Locus,
		   FILE* file,
		   exons* exons, 
		   long* nucleotides, 
		   genes* genes)
{
  char line[MAXLINE];
  char c;
  long lastPos1;
  char mess[MAXSTRING];
 
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;

  char lineTmp[MAXSTRING];

  long i;
  long iExon[STRANDS];
  int index;

  char lastGroup[MAXSTRING];
  char sGroup[MAXSTRING];

  /* Coments: line begins with # */
  /* gff format = "Name  Source  Type  Begin  End  Score  Strand  Frame group */
  i = 0;
  iExon[FORWARD] = 0; 
  iExon[REVERSE] = 0; 
  lastPos1 = -INFI;
  *nucleotides = 0;
  strcpy(lastGroup,NOGROUP);
  genes->numGenes[FORWARD] = NOGENES;
  genes->numGenes[REVERSE] = NOGENES;

  /* $ means end of current sequence of exons */
  /* Skip comment lines */
  while(fgets(line,MAXLINE,file)!=NULL && !(line[0]=='#' && line[1]=='$'))
    {
      if((line[0]=='#' && line[1]!='$') || line[0]=='\n')
		{
		  /* Skip this line */
		  sprintf(mess,"Skipping comment line");
		  printMess(mess);
		}
      else
		{
		  strcpy(lineTmp,line);

		  /* For each line extract the features (GFF format) */
          line1 = (char *) strtok(line,"\t");
          line2 = (char *) strtok(NULL,"\t");
          line3 = (char *) strtok(NULL,"\t"); 
          line4 = (char *) strtok(NULL,"\t");
		  line5 = (char *) strtok(NULL,"\t");
		  line6 = (char *) strtok(NULL,"\t");
          line7 = (char *) strtok(NULL,"\t");
          line8 = (char *) strtok(NULL,"\t");
          line9 = (char *) strtok(NULL,"\n");
		  
		  if (line1 == NULL || line2 == NULL || line3 == NULL ||
			  line4 == NULL || line5 == NULL || line6 == NULL ||
			  line7 == NULL || line8 == NULL)
			{
			  sprintf(mess, "Bad format: Exon GFF %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  
		  /* According to the strand, store in exons[+] or exons[-] */
		  if (sscanf(line7,"%c",&c)!= 1)
			{
			  sprintf(mess, "Bad format Strand: Exon %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  else
			{
			  /* select strand */
			  index = (c == '+')? FORWARD : REVERSE;
			  (exons->exon[index] + iExon[index])->Strand = c;
			}
		  
		  /* 1. Locus matching  */
		  if (strcmp(Locus,line1))
			{
			  sprintf(mess, "Mismatch between locus: %s (real) -- %s (pred)\n",
					  Locus,
					  line1);
//			  printError(mess);
//			  printMess(mess);
			}
		  
		  /* 3. Exon Type */
		  if (sscanf(line3,"%s",(exons->exon[index] + iExon[index])->Type) != 1)
			{
			  sprintf(mess, "Bad format Type: Exon %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  
		  /* 4. Left position */
		  if (sscanf(line4,"%ld",&((exons->exon[index] + iExon[index])->Position1)) != 1)
			{
			  sprintf(mess, "Bad format Position1: Exon %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  
		  /* 4.b. Filename sort by position */
		  if ((exons->exon[index] + iExon[index])->Position1 <= lastPos1)
			{
			  sprintf(mess, "Bad position(not sorted): Exon %ld (%s)\n",i,lineTmp);
			  printError(mess); 
			}
		  lastPos1 = (exons->exon[index] + iExon[index])->Position1;
		  
		  /* 5. Right position */
		  if (sscanf(line5,"%ld",&((exons->exon[index] + iExon[index])->Position2)) != 1)
			{
			  sprintf(mess, "Bad format Position2: Exon %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  
		  /* 6. Score = '.' or float */
		  if (sscanf(line6,"%lf",&((exons->exon[index] + iExon[index])->Score)) != 1)
			{
			  if ((sscanf(line6,"%c",&c)!= 1) || (c!='.'))
				{
				  sprintf(mess, "Bad format Score: Exon %ld (%s)\n",i,lineTmp);
				  printError(mess);
				}
			  (exons->exon[index] + iExon[index])->Score = NOSCORE;
			}
		  
		  /* 7. Strand was done before */
	  	  
		  /* 8. Frame = '.' or integer */
		  if (sscanf(line8,"%hd",&((exons->exon[index] + iExon[index])->Frame)) != 1)
			{
			  if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
				{
				  sprintf(mess, "Bad format Frame: Exon %ld (%s)\n",i,lineTmp);
				  printError(mess);
				}
			  (exons->exon[index] + iExon[index])->Frame = NOFRAME;
			}	  
		  
		  /* 9. Group */
		  if (line9 != NULL)
			{
			  if (sscanf(line9,"%s",sGroup) != 1)
				{
				  sprintf(mess, "Bad format Group: Exon %ld (%s)\n",i,lineTmp);
				  printError(mess);
				}
			}
		  else
			{
			  sprintf(mess, "Group not found: Exon %ld (%s)\n",i,lineTmp);
			  printError(mess);
			}
		  
		  /* Post-processing of the read exon */
		  /* a. Computing total number of nucleotides */
		  *nucleotides += (exons->exon[index] + iExon[index])->Position2 - (exons->exon[index] + iExon[index])->Position1 + 1;
		  
		  /* b. Building genes information */
		  /* The same gene (group)...? */
		  if (!strcmp(sGroup, lastGroup))
			{
			  (genes->gen[index] + genes->numGenes[index])->numExons++;
			  (genes->gen[index] + genes->numGenes[index])->Last = (exons->exon[index] + iExon[index]);
			}
		  else
			{
			  /* Init new gene */
			  genes->numGenes[index]++;
			  
			  (genes->gen[index] + genes->numGenes[index])->numExons = 1;
			  
			  (genes->gen[index] + genes->numGenes[index])->First = (exons->exon[index] + iExon[index]);
			  (genes->gen[index] + genes->numGenes[index])->Last = (exons->exon[index] + iExon[index]);
			}
		  
		  /* Updating lastGroup */
		  strcpy(lastGroup,sGroup);
		  
		  /* Process features from current exon */
		  if (i > NUMEXONS)
			printError("Too many records: Change NUMEXONS parameter");
		  
		  i++;
		  iExon[index]++;
		}
    }/* end of while */
  
  /* Updating number of genes */
  if (genes->numGenes[FORWARD] == NOGENES)
    genes->numGenes[FORWARD] = 0;
  else
    genes->numGenes[FORWARD]++;

  if (genes->numGenes[REVERSE] == NOGENES)
    genes->numGenes[REVERSE] = 0; 
  else
    genes->numGenes[REVERSE]++;
  
  /* Updating number of exons */
  for (index=0; index < STRANDS; index++)
    exons->numExons[index] = iExon[index];

  return(i);
}


