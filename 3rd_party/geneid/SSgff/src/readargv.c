/*************************************************************************
*                                                                        *
*   Module: readargv                                                     *
*                                                                        *
*   Read set up options and filenames from user input                    *
*                                                                        *
*   This file is part of the geneid 1.1 distribution                     *
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

/*  $Id: readargv.c,v 1.10 2001/12/18 16:23:47 eblanco Exp $  */

#include "SSgff.h"

/* SSgff.c external vars */
extern int 	SFP, SDP, SAP, STP, ALLS,
            EFP, EIP, ETP, ESP, ALLE, INT,
            CDS, TRN,
            VRB;

/* required by getopts */
extern char* optarg;
extern int optind;

char* USAGE="Incorrect usage:\nNAME\n\tSSgff - a program to extract features from  genomic sequences\nSYNOPSIS\n\tSSgff\t[-bdaefitxsz][-h]\n\t\t<locus_seq_in_fasta_format> <features_in_gff_format>\n\n";

void printHelp()
{
  printf("\n\tSSgff: Setup options\n");
  printf("\t------------------------\n\n");
  
  printf("\t-b: Output Start codons\n");
  printf("\t-d: Output Donor splice sites\n");
  printf("\t-a: Output Acceptor splice sites\n");
  printf("\t-e: Output Stop codons\n");

  printf("\t-s: Output Sites\n");  

  printf("\t-F: Output Initial exons\n");
  printf("\t-I: Output Internal exons\n");
  printf("\t-T: Output Terminal exons\n");
  printf("\t-S: Output Single genes\n");

  printf("\t-E: Disables default output: all exons \n");
  
  printf("\t-i: Output Introns\n");

  printf("\t-c: Output complete CDS\n");
  printf("\t-t: Output primary trancripts\n\n");
  
  
}

void readargv (int argc,char* argv[],
			   char* SequenceFile,
			   char* ExonsFile) 
{
  int c;
  int error=0;
  int printOptions =0;
  char mess[MAXSTRING];
  
  /* Reading setup options */
  while ((c = getopt(argc,argv,"bdaeFITSEictvhs")) != -1)
    switch(c)
      {
      case 'b': SFP++;
		printOptions++;
		break;
      case 'd': SDP++;
		printOptions++;
		break;
      case 'a': SAP++;
		printOptions++;
		break;
      case 'e': STP++;
		printOptions++;
		break;
	  case 's': ALLS++;
		printOptions++;
		break;
      case 'F': EFP++;
		printOptions++;
		break;
      case 'I': EIP++;
		printOptions++;
		break;
      case 'T': ETP++;
		printOptions++;
		break;
	  case 'S': ESP++;
		printOptions++;
		break;
      case 'E': ALLE++;
		printOptions++;
		break;
	  case 'i': INT++;
		printOptions++;
		break;
      case 'c': CDS++;
		printOptions++;
		break;
      case 't': TRN++;
		printOptions++;
		break;
      case 'v': VRB++;
		break;
      case '?':error++;
		break; 
      case 'h': printHelp();
		exit(1);
		break;
      }
  
  /* Setup Errors (a): Incompatible options selected */
  //  if (SAP && ESP)
  //  printError("Incompatible options( -D | -G)");
  
  if (error)
    printError(USAGE);
  
  /* Setup Errors (b): Wrong number of filenames */
  /* Read the name of the input fasta file */
  if (optind < argc)
    {
      strcpy(SequenceFile,argv[optind]);
      optind++;
    }
  else
	{
	  sprintf(mess,"Two filename is needed: DNA sequence(fasta format) and features (GFF format)\n%s",USAGE);
	  printError(mess);
	}

  if (optind < argc)
    {
      strcpy(ExonsFile,argv[optind]);
      optind++;
      if (optind < argc)
		{ 
          sprintf(mess,"Two  filenames required but more than two presented\n%s",USAGE);
		  printError(mess);
		}
    }
  else
	{
	  sprintf(mess,"Two filename is needed: DNA sequence(fasta format) and features (GFF format)\n%s",USAGE);
	  printError(mess);
	}

}
