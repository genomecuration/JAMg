/*************************************************************************
*                                                                        *
*   Module: SSgff                                                        *
*                                                                        *
*   Formatted output of geneid (GFF, default and extended)               *
*                                                                        *
*   This file is part of the SSgff  Distribution                         *
*                                                                        *
*     Copyright (C) 2001 - Genis PARRA FARRE                             *
*                          Enrique BLANCO GARCIA                         *
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

#include "SSgff.h"

/* SSgff setup flags */
int
  /* sites to print */
  SFP=0, SDP=0, SAP=0, STP=0, ALLS=0,
  /* exons to print */
  EFP=0, EIP=0, ETP=0, ESP=0, ALLE=0, INT=0,
  /* transcrip */
  CDS=0,TRN=0,
  /* Verbose flag (memory/processing information) */
  VRB=0;


int main (int argc, char *argv[])
{
  
  /* input files names */
  char FastaFile[FILENAMELENGTH], GffFile[FILENAMELENGTH];

  /* sequence */
  char Locus[MAXLINE];
  char *Sequence;

  /* Auxiliary arrays */
  char *AuxExon;
  char *AuxGene;
  char *AuxExon_tmp;

  /* exons array */
  exonGFF *exons;
  pack_gene *genes;

  /* variables */
  char mess[MAXSTRING];
  long numex;
  long numgen;
  long c, g;
  long length, lengthCDS;
  long pos1, pos2;
  long pos1_site, pos2_site, pos1_intron, pos2_intron ;
 

  /* Auxiliary pointers */ 
  exonGFF *AuxExonPtr;
  exonGFF *AuxExonPtr_tmp;

  /* dictionary variable for genes names*/
  dict *d;

  /* input files */
  FILE *fFastaFile;
  FILE *fGffFile;

    
  /* memory for the sequence */ 
  if ((Sequence = (char *) calloc(SEQLENGTH,sizeof(char))) == NULL)
	printError("Not enough memory: Sequence"); 
  if ((AuxExon = (char *) calloc(MAXAA * LENGTHCODON,sizeof(char))) == NULL)
	printError("Not enough memory: AuxExon");
  if ((AuxGene = (char *) calloc(MAXAA * LENGTHCODON,sizeof(char))) == NULL)
	printError("Not enough memory: AuxGene");
  if ((AuxExon_tmp = (char *) calloc(MAXAA * LENGTHCODON,sizeof(char))) == NULL)
	printError("Not enough memory: AuxExon_tmp");
  if ((exons = (exonGFF *) calloc(NUMEXONS,sizeof(exonGFF))) == NULL)
	printError("Not enough memory: exons"); 
  if ((genes = (pack_gene *) calloc(NUMEXONS,sizeof(pack_gene))) == NULL)
	printError("Not enough memory: genes"); 
  if ((AuxExonPtr_tmp = (exonGFF *)malloc(sizeof(exonGFF))) == NULL)
    printError("Not enough memory: AuxExonPtr_tmp");
  if ((AuxExonPtr = (exonGFF *)malloc(sizeof(exonGFF))) == NULL)
    printError("Not enough memory: AuxExonPtr");
  if ((d = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough memory: dictionary of exon types");

  /*  Read setup options */
  readargv(argc,argv,FastaFile,GffFile);
  printMess("\n\n\t\t\t** Executing SSgff 2001 gparra@imim.es **\n\n"); 
                                  
  /* Open files*/ 
  fFastaFile = OpenFile(FastaFile); 
  fGffFile = OpenFile(GffFile); 
 
  /* read  the sequence */  
  sprintf(mess,"Reading FASTA file: %s", FastaFile); printMess(mess);
  length = ReadSequence (fFastaFile, Locus, Sequence);
  sprintf(mess,"DNA sequence %s  readed (%ld bp)", Locus, length); printMess(mess);

  /* read exons */
  sprintf(mess,"Reading GFF file: %s ",GffFile); printMess(mess); 
  numex = ReadExonsGFF (fGffFile, exons, genes, d, &numgen, length);
  sprintf(mess,"%ld exons readed grouped on %ld genes ",numex, numgen); printMess(mess);
 

  /* Scanning and printting the subsequences */  
  for (g = 0; g < numgen; g++)
	{
	  sprintf(mess,"Extracting information from group:  %s , %ld exons",
			  (genes+g)->First->Group,(genes+g)->numExons+1); printMess(mess);
	  lengthCDS = 0;

	  /* Initializing firts exon before the loop */
	  AuxExonPtr = (genes+g)->First;

	  for (c = 0; c <= (genes+g)->numExons; c++) 
		{

		  /* For other exons except the first one */
		  if (c != 0)
			{
			  AuxExonPtr_tmp  = AuxExonPtr;
			  AuxExonPtr = AuxExonPtr_tmp->NextExon;
			  /* Intron setting part2 and Printing */
			  if (INT)
				{
				  pos2_intron = (AuxExonPtr->Position1)-1;
				  printHeadIntron(AuxExonPtr->Group, c, Locus);
				  printFasta(pos1_intron, pos2_intron, AuxExonPtr->Strand, Sequence);
				}
			}
		  /* Setting intron pos1 for the next loop */
		  pos1_intron = (AuxExonPtr->Position2)+1;

		  /* Defining exon position1 and position2 */
		  pos1 = AuxExonPtr->Position1;
		  pos2 = AuxExonPtr->Position2;

          CheckBoundaries(&pos1, &pos2, length);
		  if (pos2 - pos1 + 1 > MAXAA * LENGTHCODON)
			printError("Not enough memory to hold exons : change MAXAA parameter");

		  lengthCDS += pos2 - pos1 + 1;
		  if (lengthCDS > MAXAA * LENGTHCODON)
			printError("Not enough memory to hold cds : change MAXAA parameter");
	  
		  /* Printing all exons: DEFAULT */
		  if (!(ALLE))
			{
			  printHeadExon(AuxExonPtr->Group, c, Locus);
			  printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
			}
		  
		  /* Joining all the CDS fragments */
		  if (CDS)
			{
			  if ((AuxExonPtr)->Strand == '-') 
				{
				  ReverseSubSequence( pos1 - 1, pos2 - 1, Sequence, AuxExon);
				  AuxExon[ pos2 - pos1 + 1] = '\0';
				  /* Joining the exons to build the complete CDS in reverse strand*/
				  strcpy (AuxExon_tmp,AuxExon);
				  strcat(AuxExon_tmp,AuxGene);
				  strcpy (AuxGene,AuxExon_tmp);
				}
			  else
				{
				  strncpy(AuxExon, 
						  Sequence + pos1 - 1,
						  pos2 - pos1 + 1);
				  AuxExon[pos2 - pos1 + 1] = '\0';
				  /* Joining the exons to build the complete CDS in forward strand*/
				  strcat(AuxGene,AuxExon);
				} 
			}

		  /* Printing sites or exons selected by type*/

		  /* Reverse strand */
		  if (AuxExonPtr->Strand == '-') 
			{
			  /* First Exons */
			  if (!strcmp(AuxExonPtr->Type, SFIRST))
				{
				  /* Printing only first*/
				  if (EFP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* START Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos2 - DOWN_START;
					  pos2_site = pos2 + UP_START;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus, SSTART);
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				  /* DONOR Site */
				  if (ALLS || SDP)
					{
					  pos1_site = pos1 - DOWN_DONOR;
					  pos2_site = pos1 + UP_DONOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus, SDONOR);
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				}

			  /* Internal Exons */
			  else if (!strcmp(AuxExonPtr->Type, SINTERNAL))
				{
				  /* Printing only internal*/ 
				  if (EIP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* ACCEPTOR Site */
				  if (ALLS || SAP)
				  {
					pos1_site = pos2 - DOWN_ACCEPTOR;
					pos2_site = pos2 + UP_ACCEPTOR;
					CheckBoundaries(&pos1_site, &pos2_site, length);
					printHeadSite(AuxExonPtr->Group, c, Locus,  SACCEPTOR);
					printFasta(pos1_site, pos2_site, '-', Sequence);
				  }
				  /* DONOR Site */
				  if (ALLS || SDP)
					{
					  pos1_site = pos1 - DOWN_DONOR;
					  pos2_site = pos1 + UP_DONOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SDONOR); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				}
			  /* Terminal Exons */
			  else if (!strcmp(AuxExonPtr->Type, STERMINAL))
				{
				  /* Printing only terminal*/ 
				  if (ETP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}

				  /* Acceptor Site */
				  if (ALLS || SAP)
					{
					  pos1_site = pos2 - DOWN_ACCEPTOR;
					  pos2_site = pos2 + UP_ACCEPTOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SACCEPTOR); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				  /* Stop Site */
				  if (ALLS || STP)
					{
					  pos1_site = pos1 - DOWN_STOP;
					  pos2_site = pos1 + UP_STOP;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTOP); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				}

			  /* Single Exons */
			  else if (!strcmp(AuxExonPtr->Type,SSINGLE))
				{
				  /* Printing only Single*/ 
				  if (ESP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* START Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos2 - DOWN_START;
					  pos2_site = pos2 + UP_START;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTART); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				  /* STOP Site */
				  if (ALLS || STP)
					{
					  pos1_site = pos1 - DOWN_STOP;
					  pos2_site = pos1 + UP_STOP;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTOP); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				}
			  /* Not defined exon types */
			  else 	
				{
				  /* First Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos1 - DOWN_SITE2;
					  pos2_site = pos1 + UP_SITE2;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  "Site_2"); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				  /* Second Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos2 - DOWN_SITE1;
					  pos2_site = pos2 + UP_SITE1;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  "Site_1"); ;
					  printFasta(pos1_site, pos2_site, '-', Sequence);
					}
				}
			}
		  /* Forward Strand */
		  else
			{
			  /* First Exons */
			  if (!strcmp(AuxExonPtr->Type, SFIRST))
				{
				  /* Printing only first*/ 
				  if (EFP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* START Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos1 - UP_START;
					  pos2_site = pos1 + DOWN_START;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTART); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				  /* DONOR Site */
				  if (ALLS || SDP)
					{
					  pos1_site = pos2 - UP_DONOR;
					  pos2_site = pos2 + DOWN_DONOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SDONOR); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				}

			  /* Internal Exons */
			  else if (!strcmp(AuxExonPtr->Type, SINTERNAL))
				{
				  /* Printing only Internal exons*/ 
				  if (EIP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* ACCEPTOR Site */
				  if (ALLS || SAP)
					{
					  pos1_site = pos1 - UP_ACCEPTOR;
					  pos2_site = pos1 + DOWN_ACCEPTOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SACCEPTOR); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				  /* DONOR Site */
				  if (ALLS || SDP)
				  {
					pos1_site = pos2 - UP_DONOR;
					pos2_site = pos2 + DOWN_DONOR;
					CheckBoundaries(&pos1_site, &pos2_site, length);
					printHeadSite(AuxExonPtr->Group, c, Locus,  SDONOR); 
					printFasta(pos1_site, pos2_site, '+', Sequence);
				  }
				}
			  /* Terminal Exons */
			  else if (!strcmp(AuxExonPtr->Type, STERMINAL))
				{
				  /* Printing only terminal*/ 
				  if (ETP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1, pos2, (AuxExonPtr)->Strand, Sequence);
					}
				  /* ACCEPTOR Site */
				  if (ALLS || SAP)
					{
					  pos1_site = pos1 -  UP_ACCEPTOR;
					  pos2_site = pos1 + DOWN_ACCEPTOR;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SACCEPTOR); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				  /* STOP Site */
				  if (ALLS || STP)
					{
					  pos1_site = pos2 - UP_STOP ;
					  pos2_site = pos2 + DOWN_STOP;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTOP); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				}

			  /* Single Exons */
			  else if (!strcmp(AuxExonPtr->Type,SSINGLE))
				{
				  /* Printing only first*/ 
				  if (ESP)
					{
					  printHeadExon(AuxExonPtr->Group, c, Locus);
					  printFasta(pos1,pos2,(AuxExonPtr)->Strand,Sequence);
					}
				  /* START Site */
				  if (ALLS || SFP)
					{
					  pos1_site = pos1 - UP_START;
					  pos2_site = pos1 + DOWN_START;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTART); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				  /* STOP Site */
				  if (ALLS || STP)
					{
					  pos1_site = pos2 - UP_STOP;
					  pos2_site = pos2 + DOWN_STOP;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus,  SSTOP); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				}
			  else 	
				{
				  /* Print sites if Exon is not defined */
				  /* First Site */
				  if (ALLS)
					{
					  pos1_site = pos1 - UP_SITE1;
					  pos2_site = pos1 + DOWN_SITE1;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus, "Site_1"); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				  /* Second Site */
				  if (ALLS)
					{
					  pos1_site = pos2 - UP_SITE2;
					  pos2_site = pos2 + DOWN_SITE2;
					  CheckBoundaries(&pos1_site, &pos2_site, length);
					  printHeadSite(AuxExonPtr->Group, c, Locus, "Site_2"); 
					  printFasta(pos1_site, pos2_site, '+', Sequence);
					}
				}  
			}
			
		}
	  
	  /* Printing complete CDS */
	  if (CDS)
		{
		  printFastaGene(AuxExonPtr->Group, Locus, AuxGene);
		  /* Initialize AuxGene array */
		  AuxGene[0] = '\0'; 
		}
	  /* Printing Transcripts*/
	  if (TRN) 
		{
		  /* Defining exon position1 and position2 */
		  pos1 = ((genes+g)->First->Strand == '+') ? 
			(genes+g)->First->Position1 - UP_TRANSCRIP : (genes+g)->First->Position1 - DOWN_TRANSCRIP;
		  pos2 =  ((genes+g)->First->Strand == '+') ?
			(genes+g)->Last->Position2 + DOWN_TRANSCRIP : (genes+g)->Last->Position2 + UP_TRANSCRIP;
		  
		  CheckBoundaries(&pos1, &pos2, length);
		  
		  /* Extract the squence and stored them in AuxTrans */ 
		  printHeadTranscript ((genes+g)->First->Group, Locus);  
		  printFasta (pos1, pos2, (genes+g)->First->Strand, Sequence);
		}
	}
  return(0);
} 



