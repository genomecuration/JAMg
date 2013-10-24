/*************************************************************************
*                                                                        *
*   Module: readargv                                                     *
*                                                                        *
*   Get setup options and filenames of user input.                       *
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
extern int AST;
extern int TST;
extern int SST;

extern char* optarg;
extern int optind;

char *USAGE="Incorrect usage.\nNAME\n\tevaluation - a program to measure gene prediction accuracy\nSYNOPSIS\n\tevaluation [-vats] <Predict_genes> <Real_genes>";

void printHelp()
{
  printf("Short Manual for evaluation:\n");
  printf("------------------------\n\n");
  printf("Setup Options:\n\n");
 
  printf("\t-v: Verbose. Print all messages\n");
  printf("\t-a: Average. Print average stats (more than 1 sequence)\n");
  printf("\t-t: Total. Print total stats (more than 1 sequence)\n");
  printf("\t-s: Short. Print a short output\n");

  printf("\t-h: Show this Short handbook\n");
}

void readargv (int argc,char* argv[],
	       char* PredictionsFile,
	       char* RealFile) 
{
  int c;
  char mess[MAXSTRING];
  int error=0;

  /* Reading setup options */
  while ((c = getopt(argc,argv,"vsath")) != -1)
    switch(c)
      {
      case 'v': VRB++; 
	break;
      case 's': SST++; 
	break;
      case 'a': AST++; 
	break;
      case 't': TST++; 
	break;
      case '?':error++;
	break; 
      case 'h': printHelp();
	exit(1);
	break;
      }

  if (error)
    printError(USAGE);

  /* Setup Errors: Wrong number of filenames */
  /* 1. Get the name of the predictions file */
  if (optind < argc)
    {
      strcpy(PredictionsFile,argv[optind]);
      optind++;
      /* 2. Get the name of the annotations file */
      if (optind < argc)
	{
	  strcpy(RealFile,argv[optind]);
	  optind++;
	  if (optind < argc)
	    { 
	      sprintf(mess,"Too many files. Only 2 filenames needed\n%s",USAGE);
	      printError(mess);
	    }
	}
      else
	{
	  sprintf(mess,"Where is the real Annotations file?\n%s",USAGE);
	  printError(mess);
	}
    }
  else
    {
      sprintf(mess,"Where is the Predictions file?\n%s",USAGE);
      printError(mess);
    }
}










