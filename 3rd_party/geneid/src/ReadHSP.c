/*************************************************************************
*                                                                        *
*   Module: ReadHSP                                                      *
*                                                                        *
*   Loading High-score Segment Pairs (protein alignment) in GFF format   *
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

/*  $Id: ReadHSP.c,v 1.6 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* According to Locusname, select a group of HSPs to be sorted */
packHSP* SelectHSP(packExternalInformation* external,
				   char* Locus,
				   long LengthSequence)
{
  int a;
  int frame;
  long i;
  long pos1,pos2;
  packHSP* p;
  
  /* 1. Select the position in the array corresponding to Locus */
  a = getkeyDict(external->locusNames, Locus);
  if (a == NOTFOUND)
	p = NULL;
  else
	{
	  p = external->homology[a];
  
	  /* Recomputing positions in HSPs from reverse strand */
	  /* The visited field prevents from repeating the sorting */
	  if (!p->visited)
		{
		  p->visited = 1;
		  
		  /* 2. Recompute HSPs positions in reverse strand: frames 3, 4 and 5 */
		  for (frame=FRAMES; frame < 2*FRAMES; frame++)
			{
			  /* Calcule new positions */
			  for(i=0; i < p->nSegments[frame]; i++)
				{
				  pos1 = p->sPairs[frame][i]->Pos1;
				  pos2 = p->sPairs[frame][i]->Pos2;
				  
				  p->sPairs[frame][i]->Pos1 =
					LengthSequence - pos2 + 1;
				  
				  p->sPairs[frame][i]->Pos2 =
					LengthSequence - pos1 + 1;
				}
			}

		  /* 3. Quick-Sorting of FWD HSPs and RVS HSPs */
		  SortHSPs(p);
		}
	  
	  /* 4. Reset partial counter of used HSPs */
	  for(i=0; i < STRANDS*FRAMES; i++)
		external->iSegments[i] = 0;
	}

  return(p);
}

/* Input blast HSPs from external file */
long ReadHSP (char* FileName, packExternalInformation* external)
{
  long i,j;
  FILE* file;
  char line[MAXLINE];
  char lineCopy[MAXLINE];
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;

  /* Identifier for a sequence (from dictionary) */
  int a;

  /* If frame = '.' then make three copies of current SR (3 frames) */
  int three;
  char Locusname[LOCUSLENGTH];
  char mess[MAXSTRING];
  long pos1, pos2;
  float score;
  char strand;
  short frame;
  char c;
  
  
  /* Coments: line begins with # */
  /* gff format = Name  Source  Type  Begin  End  Score  Strand  Frame */
  /* group is allowed but skipped */
  
  /* 0. Opening the HSP file */
  if ((file=fopen(FileName, "r"))==NULL)
    printError("The homology file can not be opened to read");
  
  /* 1. Read HSPs */
  i = 0;
  three = 0;
  for(i=0; i<MAXNSEQUENCES; i++)
  
  while(fgets(line,MAXLINE,file)!=NULL)
    {
      /* 2.a. Comment or empty line: forget it */
      if (line[0]=='#' || line[0]=='\n')
		printMess("Skipping comment line in HSPs file");
      else
		{
          /* Backup copy of the line to show error messages */
          strcpy(lineCopy,line);
		  
          /* Extracting GFF features */
          line1 = (char *) strtok(line,"\t");
          line2 = (char *) strtok(NULL,"\t");
          line3 = (char *) strtok(NULL,"\t"); 
          line4 = (char *) strtok(NULL,"\t");
          line5 = (char *) strtok(NULL,"\t");
          line6 = (char *) strtok(NULL,"\t");
          line7 = (char *) strtok(NULL,"\t");
          line8 = (char *) strtok(NULL,"\t");
          line9 = (char *) strtok(NULL,"\n");
		  
		  /* There are 8 mandatory columns and the last one is optional */
		  if (line1 == NULL || line2 == NULL || line3 == NULL ||
			  line4 == NULL || line5 == NULL || line6 == NULL ||
			  line7 == NULL || line8 == NULL)
			{
			  sprintf(mess,"Wrong GFF format in HSPs (number of records):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 1. Locusname: leave the exon into the correct array */
		  if (sscanf(line1,"%s",Locusname) != 1)
			{
			  sprintf(mess,"Wrong GFF format in HSPs (locusname):\n-->%s\n",lineCopy);
			  printError(mess);
			}
	  
		  /* Look-up the ID for that sequence */
		  a = setkeyDict(external->locusNames,Locusname);
		  if (a >= MAXNSEQUENCES)
			printError("Too many DNA sequences: increase MAXNSEQUENCES parameter");
		  
		  /* 4. Starting position */
		  if (sscanf(line4,"%ld",&pos1) != 1)
			{
			  sprintf(mess,"Wrong GFF format in HSPs (starting position):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 5. Finishing position */
		  if (sscanf(line5,"%ld",&pos2) != 1)
			{
			  sprintf(mess,"Wrong GFF format in HSPs (finishing position):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 6. HSP score */
		  if (sscanf(line6,"%f",&score) != 1)
			{
			  sprintf(mess,"Wrong GFF format in HSPs (score):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 7. Strand (reading sense) [+|-] */
		  if ((sscanf(line7,"%c",&strand)!= 1) ||
			  ((strand != '+') && (strand != '-')))
			{       
			  sprintf(mess,"Wrong GFF format in HSPs (strand):\n-->%s\n",lineCopy);
			  printError(mess);
			}
		  
		  /* 8. Frame = integer or '.' */
		  if (sscanf(line8,"%hd",&frame) != 1)
			{
			  /* Is it a dot? */
			  if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
				{
				  sprintf(mess,"Wrong GFF format in HSPs (frame):\n-->%s\n",lineCopy);
				  printError(mess);
				}
			  /* make three copies for this SR */
			  three = 1;
			}
		  else
			{
			  /* Checking input frame between 1..3 */
			  if ((frame < 1) || (frame > 3))
				{
				  sprintf(mess,"Wrong GFF value in HSPs frame (must be in [1..3]):\n-->%s\n",
						  lineCopy);
				  printError(mess);
				}
			}

		  /* Allocating current hsp in packHSP */
		  /* a) Forward sense */
		  if (strand == '+')
			{
			  /* a1. Unknown frame: 3 copies */
			  if (three)
				{
				  /* Replication: 3 SRs for the same coordinates */
				  for (frame=0; frame < FRAMES; frame++)
					{
					  /* New item HSP */
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]
						= (HSP*) RequestNewHSP();
					  
					  /* Setting SR data */
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos1 = pos1;
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos2 = pos2;
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Score = score;
					  external->homology[a]->nSegments[frame]++;
					}
				}
			  /* a2. Known frame */
			  else
				{
				  /* Convert blast frame 3 into frame 0 to store it */
				  frame = (frame % 3);

				  /* New item HSP */
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]
					= (HSP*) RequestNewHSP();
				  
				  /* Setting HSP data */
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos1 = pos1;
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos2 = pos2;
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Score = score;
				  
				  external->homology[a]->nSegments[frame]++;
				}
			} /* end FWD */
		  /* b) Reverse sense */
		  else
			{
			  /* b1. Unknown frame: 3 copies */
			  if (three)
				{
				  /* Replication: 3 HSPs for the same coordinates */
				  /* Frames 3, 4 and 5 as in 0, 1 and 2 (reverse) */
				  for (frame=FRAMES; frame < 2*FRAMES; frame++)
					{
					  /* New item HSP */
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]
						= (HSP*) RequestNewHSP();

					  /* Setting HSP data */
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos1 =
                        pos1;
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos2 =
                        pos2;
					  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Score = 
						score;
					  external->homology[a]->nSegments[frame]++;
					}
				}
			  /* b2. Known frame */
			  else
				{
				  /* Convert blast frame 3 into frame 0 to store it */
				  /* Frames 3, 4 and 5 equal 0, 1 and 2 (reverse) */
				  frame = (frame % 3);
				  frame = FRAMES + frame;

				  /* New item HSP */
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]
					= (HSP*) RequestNewHSP();
				  
				  /* Setting HSP data */
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos1 =
					pos1;
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Pos2 =
					pos2;
				  external->homology[a]->sPairs[frame][external->homology[a]->nSegments[frame]]->Score = 
					score;
				  external->homology[a]->nSegments[frame]++;
				}
			} /* end RVS */
  
		  /* ready for next HSP... */
		  external->homology[a]->nTotalSegments = (three)? 
			external->homology[a]->nTotalSegments + 3: 
			external->homology[a]->nTotalSegments + 1;
		  three = 0;

		  if ((external->homology[a]->nTotalSegments + FRAMES) >= MAXHSP)
			printError("Too many HSPs: increase MAXHSP definition");

		} /* end of if-comment */
    } /* end of while*/

  fclose(file);
  
  /* Obtain the number of different sequences */
  external->nSequences = external->locusNames->nextFree;

  /* Return the number of created HSPs (including replications) */
  for(i=0,j=0; j < external->nSequences; j++) 
	i = i + external->homology[j]->nTotalSegments;


  return(i);
}
