/*************************************************************************
*                                                                        *
*   Module: readparam                                                    *
*                                                                        *
*   Reading statistic parameters and gene construction model             *
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

/*  $Id: readparam.c,v 1.17 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

extern float NO_SCORE;
extern int PPT;
extern int BP;
extern int U12GTAG;
extern int U12ATAC;
extern int U2GCAG;
extern int U2GTA;
extern int U2GTG;
extern int U2GTY;
extern int SGE;
extern int PAS;
extern int BKGD_SUBTRACT_FLANK_LENGTH;
extern short SPLICECLASSES;
extern float U12_SPLICE_SCORE_THRESH;
extern float U12_EXON_SCORE_THRESH;
extern float U12EW;
extern float RSSMARKOVSCORE;
extern float RSSDON;
extern float RSSACC;
extern float EvidenceEW;
extern float EvidenceFactor;

/* Numeric values: read one line skipping comments and empty lines */
void readLine(FILE *File, char* line)
{
  char* res;
 
  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#')|| (line[0]=='\n')))
    res = fgets(line,MAXLINE,File);

  if (res == NULL)
    printError("Parameter file: unexpected end of reading");
}

/* Headers of values: read one line skipping comments and empty lines */
void readHeader(FILE* File, char* line)
{
  char* res;

  res = fgets(line,MAXLINE,File);
  while((res !=NULL) && ((line[0]=='#') || (line[0]=='\n')))
	res = fgets(line,MAXLINE,File);
  
  if (res == NULL)
    printError("Parameter file: unexpected end of reading");
}

void SetProfile(profile* p, FILE* RootFile, char* signal)
{
  int i,j,x,y;
  char line[MAXLINE];
  char mess[MAXSTRING];

  /* According to the order of Markov chain, select a different method */
  /* Position weight array: transition probabilities in every position */
  switch(p->order)
	{
	case 0:
	  /* 5 combinations / pos */
	  for(i=0; i < p->dimension; i++)  
		{
		  j = 0;

		  /* Reading A,C,G and T: creating N transition */
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",
					  signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  readLine(RootFile,line);
		  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3])))!=1)
			{
			  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
			  printError(mess);
			}
		  /* Generating the N value */
		  p->transitionValues[i][j+4] = 
			(p->transitionValues[i][j] +
			p->transitionValues[i][j+1] +
			p->transitionValues[i][j+2] +
			p->transitionValues[i][j+3]) / 4;
		}
	  break;
	case 1:
	  /* 25 combinations / pos */
	  for(i=0; i < p->dimension; i++)  
		{
		  /* Reading AX,CX,GX and TX: creating XN */
		  for(j=0; j < p->dimensionTrans-5; j=j+5)
			{
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  readLine(RootFile,line);
			  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3])))!=1)
				{
				  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
				  printError(mess);
				}
			  /* Generating the XN value every iteration: AN,CN,GN,TN */
			  p->transitionValues[i][j+4] = 
				(p->transitionValues[i][j] +
				 p->transitionValues[i][j+1] +
				 p->transitionValues[i][j+2] +
				 p->transitionValues[i][j+3]) / 4;
			}
		  
		  /* Creating NA,NC,NG,NT transition values */
		  for(x=0; x < 4; x++)
			{
			   p->transitionValues[i][j+x] =  
				 (p->transitionValues[i][x] +
				  p->transitionValues[i][5+x] +
				  p->transitionValues[i][10+x] +
				  p->transitionValues[i][15+x]) / 4;
			}

		  /* Creating the value NN */
		  p->transitionValues[i][j+4] =
			(p->transitionValues[i][j] +
			 p->transitionValues[i][j+1] +
			 p->transitionValues[i][j+2] +
			 p->transitionValues[i][j+3]) / 4;
		}
	  break;
	case 2:
	  /* 125 combinations / pos: there are "dimension" positions */
	  for(i=0; i < p->dimension; i++) 
		{
		  /* Reading AXX,CXX,GXX,TXX and creating ANX,CNX,GNX,TNX */
		  for(j=0; j < p->dimensionTrans-25; j=j+25)
			{
			  /* 20 = 5^order - 5 */
			  for(y=0; y < 20; y=y+5)
				{
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+1+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+2+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  readLine(RootFile,line);
				  if ((sscanf(line,"%*d %*s %f", &(p->transitionValues[i][j+3+y])))!=1)
					{
					  sprintf(mess,"Wrong format: Transition values in %s profile",signal);
					  printError(mess);
					}
				  /* Generating the XXN value (5th value): XAN,XCN,XGN,XTN */
				  p->transitionValues[i][j+4+y] = 
					(p->transitionValues[i][j+y] +
					 p->transitionValues[i][j+1+y] +
					 p->transitionValues[i][j+2+y] +
					 p->transitionValues[i][j+3+y]) / 4;
				}
			  
			  /* Creating XNA,XNC,XNG,XNT */
			  for(x=0; x < 4; x++)
				{
				  p->transitionValues[i][j+y+x] =  
					(p->transitionValues[i][j+x] +
					 p->transitionValues[i][j+5+x] +
					 p->transitionValues[i][j+10+x] +
					 p->transitionValues[i][j+15+x]) / 4;
				}
			  
			  /* Creating the value XNN */
			  p->transitionValues[i][j+4+y] =
				(p->transitionValues[i][j+y] +
				 p->transitionValues[i][j+1+y] +
				 p->transitionValues[i][j+2+y] +
				 p->transitionValues[i][j+3+y]) / 4;			  
			}
		  
		  /* Creating NAX,NCX,NGX,NTX (j=100)*/
		  for(y=0; y < 20; y=y+5)
			{
			  /* Creating NAX... */
			  for(x=0; x < 4; x++)
				{
				  p->transitionValues[i][j+y+x] =  
					(p->transitionValues[i][y+x] +
					 p->transitionValues[i][25+y+x] +
					 p->transitionValues[i][50+y+x] +
					 p->transitionValues[i][75+y+x]) / 4;
				}
			  /* Computing NAN (x=4)*/
			  p->transitionValues[i][j+y+x] = 
				(p->transitionValues[i][j+y] + 
				 p->transitionValues[i][j+y+1] + 
				 p->transitionValues[i][j+y+2] +
				 p->transitionValues[i][j+y+3]) / 4;
			}

		  /* Creating NNX (j=100,y=20)*/
		  for(x=0; x < 4; x++)
			{
			  p->transitionValues[i][j+y+x] = 
				(p->transitionValues[i][j+x] +
				 p->transitionValues[i][j+5+x] +
				 p->transitionValues[i][j+10+x] +
				 p->transitionValues[i][j+15+x]) /4;
			}

		  /* Finally, creating NNN (j=100,y=20,x=4) */
		  p->transitionValues[i][j+y+x] = 
			(p->transitionValues[i][j+y] +
			 p->transitionValues[i][j+y+1] +
			 p->transitionValues[i][j+y+2] +
			 p->transitionValues[i][j+y+3]) /4;
		}
	  break;
	default:
	  printError("Sorry, Markov order higher than 2 not implemented yet");
	  break;
	}
}

/* Read information useful to predict Start codons, donors, stop codons */
/* H parameter: 0 for not reading header, 1 for reading the header */
void ReadProfile(FILE* RootFile, profile* p, char* signal, int H)
{
  char line[MAXLINE];
  char mess[MAXSTRING];
  int numProfileParams = 0;
  p->afactor = 0;
  p->bfactor = 1;
  p->acc_context = ACCEPTOR_CONTEXT;
  p->dist = MIN_U12BPACC_DIST;
  p->opt_dist = OPT_U12BP_DIST;
  p->penalty_factor = (float) U12BP_PENALTY_SCALING_FACTOR;

  /* Definition parameters: Length, offset, cutoff and order (Markov chain) */
  if (H==1)
	readHeader(RootFile,line);

  readLine(RootFile,line);
  numProfileParams = sscanf(line,"%d %d %f %d %f %f %d %d %d %f", 
			    &(p->dimension),
			    &(p->offset), 
			    &(p->cutoff),
			    &(p->order),
			    &(p->afactor),
			    &(p->bfactor),
			    &(p->acc_context),
			    &(p->dist),
			    &(p->opt_dist),
			    &(p->penalty_factor)
			    );
  if (numProfileParams<4)
	{
	  sprintf(mess,"Wrong format: Definition parameters in %s profile",signal);
	  printError(mess);
	}
  if (numProfileParams<6)
	{
	  p->afactor = 0;
	  p->bfactor = 1;
	}
   if (numProfileParams<6)
	{
	  p->afactor = 0;
	  p->bfactor = 1;
	}
   if (numProfileParams<6)
	{
	  p->afactor = 0;
	  p->bfactor = 1;
	}

  /* Memory to allocate the data with these fixed dimensions */
  RequestMemoryProfile(p);
  
  /* Useful to check everything is OK */
  PrintProfile(p,signal);
  
  /* Prepare and read values */
  SetProfile(p,RootFile,signal);
}

/* Read information useful to predict Acceptor splice sites */
void ReadProfileSpliceSites(FILE* RootFile, gparam* gp)
{
  char line[MAXLINE];
  char mess[MAXSTRING];
  char header[MAXSTRING];
  profile* p;
  int u12bp=0;
  int u12gtagAcc=0;
  int u12atacAcc=0;
  int u12gtagDon=0;
  int u12atacDon=0;


  /* A. Optional profiles: U12GTAG, U12ATAC Donor, acceptor and branch points 
  and U2 branch points and Poly Pyrimidine Tract */
  
  readHeader(RootFile,line);
  if ((sscanf(line,"%s",header))!=1)
	{
	  sprintf(mess,"Wrong format: header in optional profile for splice sites");
	  printError(mess);
	}

  while(strcasecmp(header,sprofileACC))
	{ 
	
	/* Read optional profiles: sprofilePPT,sprofileBP,sprofileU12BP,sprofileU12gtagACC,sprofileU12atacACC */ 

	  if (!strcasecmp(header,sprofileU12BP))
		{
		  	u12bp++;

		  	/* Reading the U12 Branch Point profile */
		  	p = gp->U12BranchPointProfile;
		  	ReadProfile(RootFile, p, sBP, 0);
		}
	  	else
		  {
		  	if (!strcasecmp(header,sprofileU12gtagACC))
			  {
			  	u12gtagAcc++;

			  	/* Reading the U12gtag acceptor profile */
			  	p = gp->U12gtagAcceptorProfile;
			  	ReadProfile(RootFile, p, sACC, 0);
			  }
			  else
				{
			  	  if (!strcasecmp(header,sprofileU12atacACC))
					{	
				  	  u12atacAcc++;

				  	  /* Reading the U12atac acceptor profile */
				  	  p = gp->U12atacAcceptorProfile;
				  	  ReadProfile(RootFile, p, sACC, 0);
					}
					else
					  {
						if (!strcasecmp(header,sprofilePPT))
						  {
					  		/* Switch on the acceptor prediction using PPTs */
					  		PPT++;

					  		/* Reading the Poly Pyrimidine Tract profile */
					  		p = gp->PolyPTractProfile;
					  		ReadProfile(RootFile, p, sPPT, 0);
						  }
				  		  else
							{
					   		  if (!strcasecmp(header,sprofileBP))
								{
						  		  /* Switch on the acceptor prediction using BPs */
						  		  BP++;

						  		  /* Reading the Branch Point profile */
						  		  p = gp->BranchPointProfile;
						  		  ReadProfile(RootFile, p, sBP, 0);
								}
					  			else
								  {
						  			sprintf(mess,"Wrong format: profile name %s \n\tis not admitted for acceptors [only %s, %s, %s, %s, %s or %s]",
								  	header, sprofileACC, sprofilePPT, sprofileBP, sprofileU12BP, sprofileU12gtagACC, sprofileU12atacACC);
						  			printError(mess);
								  }
							}
					  }
				}
		  }

	  /* Next profile for acceptor site */
	  readHeader(RootFile,line);
	  if ((sscanf(line,"%s",header))!=1)
		{
		  sprintf(mess,"Wrong format: header in optional profile for acceptor sites");
		  printError(mess);
		}
	}

  /* Reading the Acceptor site profile */
  p = gp->AcceptorProfile;
  ReadProfile(RootFile, p, sACC, 0);
  
  readHeader(RootFile,line);
  if ((sscanf(line,"%s",header))!=1)
    {
      sprintf(mess,"Wrong format: header in optional profile for splice sites");
      printError(mess);
    }

  /* Read optional profiles: sprofileU12gtagDON,sprofileU12atacDON */
  
  while(strcasecmp(header,sprofileDON))
    { 

      if (!strcasecmp(header,sprofileU12gtagDON))
	{

	  u12gtagDon++;

	  /* Reading the U12gtag donor profile */
	  p = gp->U12gtagDonorProfile;
	  ReadProfile(RootFile, p, sDON, 0);
	}
      else
	{
	  if (!strcasecmp(header,sprofileU12atacDON))
	    {

	      u12atacDon++;

	      /* Reading the U12atac donor profile */
	      p = gp->U12atacDonorProfile;
	      ReadProfile(RootFile, p, sDON, 0);
	    }
	  else
	    {
	      if (!strcasecmp(header,sprofileU2gcagDON))
		{

		  U2GCAG++;

		  /* Reading the U2gcag donor profile */
		  p = gp->U2gcagDonorProfile;
		  ReadProfile(RootFile, p, sDON, 0);
		}
	      else
		{
		  if (!strcasecmp(header,sprofileU2gtaDON))
		    {

		      U2GTA++;

		      /* Reading the U2gta donor profile */
		      p = gp->U2gtaDonorProfile;
		      ReadProfile(RootFile, p, sDON, 0);
		    }
		  else
		    {
		      if (!strcasecmp(header,sprofileU2gtgDON))
			{

			  U2GTG++;

			  /* Reading the U2gtg donor profile */
			  p = gp->U2gtgDonorProfile;
			  ReadProfile(RootFile, p, sDON, 0);
			}
		      else
			{
		  	  if (!strcasecmp(header,sprofileU2gtyDON))
			    {

			      U2GTY++;

			      /* Reading the U2gty donor profile */
			      p = gp->U2gtyDonorProfile;
			      ReadProfile(RootFile, p, sDON, 0);
			    }
			  else
			    {		  
			      sprintf(mess,"Wrong format: profile name %s \n\tis not admitted for donors [only %s, %s, %s, %s, %s, %s or %s]",
				      header, sprofileDON, sprofileU12gtagDON, sprofileU12atacDON,sprofileU2gcagDON,sprofileU2gtaDON,sprofileU2gtgDON,sprofileU2gtyDON);
			      printError(mess);		
			    }
	  		}
		    }
		}
	    }
	}
      /* Next profile for donor site */
      readHeader(RootFile,line);
      if ((sscanf(line,"%s",header))!=1)
	{
	  sprintf(mess,"Wrong format: header in optional profile for donor sites");
	  printError(mess);
	}
    }

  /* Reading the Donor site profile */
  p = gp->DonorProfile;
  ReadProfile(RootFile, p, sDON, 0);
  
  /* Switch on the site prediction of U12gtag and U12atac introns */
  if (u12bp && u12gtagAcc && u12gtagDon){ U12GTAG++;SPLICECLASSES++;}
  if (u12bp && u12atacAcc && u12atacDon){ U12ATAC++;SPLICECLASSES++;}
  
}

/* Read information about signal and exon prediction in one isochore */
/* - isochores are specific DNA regions according to the G+C content - */
void ReadIsochore(FILE* RootFile, gparam* gp)
{
  float lscore;       
  int OligoLength_1;
  int i,j,f;
  char line[MAXLINE];
  char mess[MAXSTRING];
  char header[MAXSTRING];
  
  /* 1. read boundaries of isochores */
  readHeader(RootFile,line);
  readLine(RootFile,line);

  if ((sscanf(line,"%d %d\n", 
			  &(gp->leftValue),
			  &(gp->rightValue)))!=2)
    printError("Wrong format: isochore boundaries (G+C percent)");

  sprintf(mess,"Isochores boundaries(min/max percentage): %d,%d", 
		  gp->leftValue,
		  gp->rightValue);
  printMess(mess); 

  /* 2. read cutoff (final score) to accept one predicted exon */
  readHeader(RootFile,line);
  readLine(RootFile,line);

  if ((sscanf(line,"%f %f %f %f %f\n",
			  &(gp->Initial->ExonCutoff),
			  &(gp->Internal->ExonCutoff),
			  &(gp->Terminal->ExonCutoff),
	                  &(gp->Single->ExonCutoff),
			  &(gp->utr->ExonCutoff)))< 4)
    printError("Wrong format: exon score cutoffs (number/type)");  

  sprintf(mess,"Exon cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f\t%9.3f",
		  gp->Initial->ExonCutoff,
		  gp->Internal->ExonCutoff,
		  gp->Terminal->ExonCutoff,
	          gp->Single->ExonCutoff,
	          gp->utr->ExonCutoff);
  printMess(mess); 

  /* 3. read cutoff (potential coding score) to accept one predicted exon */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",     
			  &(gp->Initial->OligoCutoff),
			  &(gp->Internal->OligoCutoff),
			  &(gp->Terminal->OligoCutoff),
			  &(gp->Single->OligoCutoff)))!=4)
    printError("Wrong format: potential coding score cutoffs (number/type)");  
  
  sprintf(mess,"Oligo cutoffs: \t%9.3f\t%9.3f\t%9.3f\t%9.3f",
		  gp->Initial->OligoCutoff,
		  gp->Internal->OligoCutoff,
		  gp->Terminal->OligoCutoff,
		  gp->Single->OligoCutoff);
  printMess(mess); 

  /* 4. Weight of signals in final exon score */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f %f\n",
			  &(gp->Initial->siteFactor),
			  &(gp->Internal->siteFactor),
			  &(gp->Terminal->siteFactor),
			  &(gp->Single->siteFactor),
			  &(gp->utr->siteFactor)))<4)
	printError("Wrong format: weight of signal scores (number/type)");  
  
  sprintf(mess,"Site factors: \t%9.2f\t%9.2f\t%9.2f\t%9.2f\t%9.2f",
		  gp->Initial->siteFactor,
		  gp->Internal->siteFactor,
		  gp->Terminal->siteFactor,
	          gp->Single->siteFactor,
	          gp->utr->siteFactor);
  printMess(mess); 

  /* 5. Weight of coding potential in final exon score */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f\n",
			  &(gp->Initial->exonFactor),
			  &(gp->Internal->exonFactor),
			  &(gp->Terminal->exonFactor),
			  &(gp->Single->exonFactor)))<4)
	printError("Wrong format: weight of coding potential scores (number/type)");  
  
  sprintf(mess,"Exon factors: \t%9.2f\t%9.2f\t%9.2f\t%9.2f",
		  gp->Initial->exonFactor,
		  gp->Internal->exonFactor,
		  gp->Terminal->exonFactor,
		  gp->Single->exonFactor);
  printMess(mess); 

  /* 6. Weight of homology information in final exon score */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f %f\n",
			  &(gp->Initial->HSPFactor),
			  &(gp->Internal->HSPFactor),
			  &(gp->Terminal->HSPFactor),
			  &(gp->Single->HSPFactor),
			  &(gp->utr->HSPFactor)))<4)

	printError("Wrong format: weight of homology scores (number/type)");  
  
  sprintf(mess,"HSP factors: \t\t%9.2f\t%9.2f\t%9.2f\t%9.2f\t%9.2f",
		  gp->Initial->HSPFactor,
		  gp->Internal->HSPFactor,
		  gp->Terminal->HSPFactor,
		  gp->Single->HSPFactor,
		  gp->utr->HSPFactor);
  printMess(mess); 

  /* 7. read weights to correct the score of exons after general cutoff */
  readHeader(RootFile,line);
   if ((sscanf(line,"%s",header))!=1)
		{
		  sprintf(mess,"Wrong format: header for exon weights and optional U12 score threshold");
		  printError(mess);
		}
  while(strcasecmp(header,sExon_weights)&& strcmp(header,"Exon_weigths"))
  { 
	 /* 1. Read RSSMARKOVSCORE for markov score to assign non-exonic recursively spliced elements */
	if(!strcasecmp(header,sRSSMARKOVSCORE)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(RSSMARKOVSCORE)))!=1)
			printError("Wrong format: RSSMARKOVSCORE value scores (number/type)");  

		  sprintf(mess,"RSSMARKOVSCORE: \t%9.2f",RSSMARKOVSCORE);
		  printMess(mess);
	}
	 /* 1. Read Evidence Exon Weight */
	if(!strcasecmp(header,sEVIDENCEW)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(EvidenceEW)))!=1)
			printError("Wrong format: EvidenceExonWeight value (number/type)");  

		  sprintf(mess,"EvidenceExonWeight: \t%9.2f",EvidenceEW);
		  printMess(mess);
	}
	 /* 1. Read Evidence Exon Factor */
	if(!strcasecmp(header,sEVIDENCEF)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(EvidenceFactor)))!=1)
			printError("Wrong format: EvidenceFactor value (number/type)");  

		  sprintf(mess,"EvidenceFactor: \t%9.2f",EvidenceFactor);
		  printMess(mess);
	}
	 /* 1. Read RSS_Donor_Score_Cutoff */
	if(!strcasecmp(header,sRSS_DONOR_SCORE_CUTOFF)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(RSSDON)))!=1)
			printError("Wrong format: RSSDON value scores (number/type)");  

		  sprintf(mess,"RSSDON: \t%9.2f",RSSDON);
		  printMess(mess);
	}
	 /* 1. Read RSSMARKOVSCORE for markov score to assign non-exonic recursively spliced elements */
	if(!strcasecmp(header,sRSS_ACCEPTOR_SCORE_CUTOFF)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(RSSACC)))!=1)
			printError("Wrong format: RSSACC value scores (number/type)");  

		  sprintf(mess,"RSSACC: \t%9.2f",RSSACC);
		  printMess(mess);
	}
	 /* 1. Read U12_SPLICE_SCORE_THRESH for sum of U12 donor and acceptor splice scores */
	if(!strcasecmp(header,sU12_SPLICE_SCORE_THRESH)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(U12_SPLICE_SCORE_THRESH)))!=1)
			printError("Wrong format: U12_SPLICE_SCORE_THRESH value scores (number/type)");  

		  sprintf(mess,"U12_SPLICE_SCORE_THRESH: \t%9.2f",
				  U12_SPLICE_SCORE_THRESH);
		  printMess(mess);
	}
	 /* 1. Read U12_EXON_SCORE_THRESH for sum of U12 donor and acceptor exon scores */
	if(!strcasecmp(header,sU12_EXON_SCORE_THRESH)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(U12_EXON_SCORE_THRESH)))!=1)
			printError("Wrong format: U12_EXON_SCORE_THRESH value scores (number/type)");  

		  sprintf(mess,"U12_EXON_SCORE_THRESH: \t%9.2f",
				  U12_EXON_SCORE_THRESH);
		  printMess(mess);
	}
	 /* 1. Read U12_EXON_WEIGHT, an additional exon weight that applies to exons flanking U12 introns */
	if(!strcasecmp(header,sU12_EXON_WEIGHT)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%f\n", &(U12EW)))!=1)
			printError("Wrong format: U12_EXON_WEIGHT value score (number/type)");  

		  sprintf(mess,"U12_EXON_WEIGHT: \t%9.2f",
				  U12EW);
		  printMess(mess);
	}
	readHeader(RootFile,line);
	if ((sscanf(line,"%s",header))!=1)
		{
		  sprintf(mess,"Wrong format: header for exon weights");
		  printError(mess);
		}
  }
  readLine(RootFile,line);
  if ((sscanf(line,"%f %f %f %f %f\n",
			  &(gp->Initial->ExonWeight),
			  &(gp->Internal->ExonWeight),
			  &(gp->Terminal->ExonWeight),
			  &(gp->Single->ExonWeight),
			  &(gp->utr->ExonWeight)))<4)
    printError("Wrong format: exon weight values (number/type)");  

  sprintf(mess,"Exon weights: \t%9.3f\t%9.3f\t%9.3f\t%9.3f\t%9.3f",
		  gp->Initial->ExonWeight,
		  gp->Internal->ExonWeight,
		  gp->Terminal->ExonWeight,
		  gp->Single->ExonWeight,
		  gp->utr->ExonWeight);
  printMess(mess);

 
  /* 8. Read splice site profiles */
  /* (a).start codon profile */
  ReadProfile(RootFile, gp->StartProfile , sSTA,1);

  /* (b).acceptor and donor site profiles */
  ReadProfileSpliceSites(RootFile, gp);

  /* (c).donor site profile */
  /* ReadProfile(RootFile, gp->DonorProfile , sDON,1); */

  /* (d).stop codon profile */
  ReadProfile(RootFile, gp->StopProfile , sSTO,1);
  
  /* 9. read coding potential log-likelihood values (Markov chains) */
  readHeader(RootFile,line);   
  if ((sscanf(line,"%s",header))!=1)
    {
      sprintf(mess,"Wrong format: header ");
      printError(mess);
      printMess(header);
    }


  while(strcasecmp(header,sMarkov)&&strcasecmp(header,"Markov_oligo_logs_file"))
    {
      /* printMess(header); */
      
      
      if (!strcasecmp(header,sprofilePolyA))
	{

	  PAS++;
	  printMess("Reading PolyA Signal Profile");
	  /* Reading the U2gty donor profile */
	  ReadProfile(RootFile, gp->PolyASignalProfile, sPOL, 0);
	  
	}
      
      readHeader(RootFile,line);
      if ((sscanf(line,"%s",header))!=1)
	{
	  sprintf(mess,"Wrong format: header ");
	  printError(mess);
	  printMess(header);
	}
  }      
  /* Next profile for Markov order */
  readLine(RootFile,line);

  if ((sscanf(line,"%d", &(gp->OligoLength)))!=1)
    printError("Wrong format: oligonucleotide length");

  sprintf(mess,"Oligonucleotide (word) length: %d",gp->OligoLength);
  printMess(mess);

  /* (a). Initial probability matrix */
  printMess("Reading Markov Initial likelihood matrix");

  /* Computing the right number of initial values to read */      
  gp->OligoDim=(int)pow((float)4,(float)gp->OligoLength);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
        {
          sprintf(mess, "Wrong format/nunber (%s): Initial Markov value", line);
          printError(mess);
        }
      gp->OligoLogsIni[f][i]=lscore;
    }

  /* (b). Transition probability matrix */
  printMess("Reading Markov Transition likelihood matrix");

  OligoLength_1= gp->OligoLength + 1;
  gp->OligoDim_1=(int)pow((float)4,(float)OligoLength_1);

  sprintf(mess,"Used oligo array size: %ld",gp->OligoDim_1 * 3);
  printMess(mess);

  readHeader(RootFile,line);
  for(j = 0; j < gp->OligoDim_1 * 3; j++)
    { 
      readLine(RootFile,line); 
      if ((sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3)
        {
          sprintf(mess, "Wrong format/number (%s): Transition Markov value", line);
          printError(mess);
        }
      gp->OligoLogsTran[f][i]=lscore;
    }

  /* 10. read maximum number of donors per acceptor site (BuildExons) */
  readHeader(RootFile,line);
  readLine(RootFile,line); 
  
  if ((sscanf(line,"%d", &(gp->MaxDonors)))!=1)
    printError("Bad format: MaxDonors value");   
    
  sprintf(mess,"Maximum donors by acceptor = %d\n", gp->MaxDonors);
  printMess(mess);
}

/* Read the input of statistics data model */
int readparam (char* name, gparam** isochores)
{
  FILE* RootFile;
  char* Geneid;
  char ExternalFileName[FILENAMELENGTH];

  int i;
  char line[MAXLINE];
  char mess[MAXSTRING];
  char header[MAXSTRING];
  int nIsochores;

  /* 0. Select parameters filename for reading it */
  /* Filename must be: option P, env.var GENEID or default (none) */
  Geneid=getenv("GENEID");

  /* (a) Using -P option */
  if (strcmp(PARAMETERFILE,name))
	{
	  sprintf(mess,"Loading parameter file by using -P option: %s",name);
	  if ((RootFile = fopen(name,"rb"))==NULL)
		printError("Parameter file (-P) can not be open to read");
	}
  /* (b) Using GENEID environment var */
  else
    if (Geneid)
      {
        sprintf(mess,"Loading parameter file by using GENEID (env. var): %s",Geneid);
		sprintf(ExternalFileName,"%s", Geneid);
		if ((RootFile = fopen(ExternalFileName,"rb"))==NULL) 
          printError("Parameter file (GENEID env.var) can not be open to read");
      }
  /* (c) Using default parameter file */
    else
      {
        sprintf(mess,"Loading parameter file default");
		if ((RootFile = fopen(name,"rb"))==NULL)
          printError("Parameter file (default) can not be open to read");
      }

  /* rootfile will be the parameter file handle descriptor */
  printMess(mess);

  /* 1. Read NO_SCORE penalty for nucleotides not supported by homology */
  readHeader(RootFile,line);
  readLine(RootFile,line);
  if ((sscanf(line,"%f\n",
			  &(NO_SCORE)))!=1)
	printError("Wrong format: NO_SCORE value scores (number/type)");  
  
  sprintf(mess,"NO_SCORE: \t%9.2f",
		  NO_SCORE);
  printMess(mess);
  /* BKGD_SUBTRACT_FLANK_LENGTH */


  /* 2. Read the number of isochores */
  readHeader(RootFile,line);
  if ((sscanf(line,"%s",header))!=1)
		{
		  sprintf(mess,"Wrong format: header for number of isochores");
		  printError(mess);
		}
  while(strcasecmp(header,sNUMISO))
  { 
    if(!strcasecmp(header,sBKGD_SUBTRACT_FLANK_LENGTH)){
		  readLine(RootFile,line);
		  if ((sscanf(line,"%d\n", &(BKGD_SUBTRACT_FLANK_LENGTH)))!=1)
			printError("Wrong format: BKGD_SUBTRACT_FLANK_LENGTH value (integer)");  

		  sprintf(mess,"BKGD_SUBTRACT_FLANK_LENGTH: \t%d",
				  BKGD_SUBTRACT_FLANK_LENGTH);
		  printMess(mess);
	}
	readHeader(RootFile,line);
	if ((sscanf(line,"%s",header))!=1)
		{
		  sprintf(mess,"Wrong format: header for number of isochores");
		  printError(mess);
		}
  }
  readLine(RootFile,line);
  if ((sscanf(line,"%d\n", &(nIsochores)))!=1)
	printError("Wrong format: Number of isochores");
    
  sprintf(mess,"Number of isochores: %d", nIsochores);
  printMess(mess);

  if (nIsochores > MAXISOCHORES || !(nIsochores>0))
    printError("Wrong value: number of isochores(MAXISOCHORES)");

  /* 3. Reading every one of the isochores */
  for(i=0; i<nIsochores; i++)
    {
      sprintf(mess,"Reading isochore %d", i+1);
      printMess(mess);
      ReadIsochore(RootFile,isochores[i]);
    }

  /* 4. Reading the GeneModel */
  readHeader(RootFile,line);
  
  /* Ready to update dictionary of exon types */
  resetDict(isochores[0]->D);
  printMess("Dictionary ready to acquire information");
   
  if (SGE)
	{
	  printMess("Using an internal Gene Model");
	  isochores[0]->nclass = ForceGeneModel(isochores[0]->D,
											isochores[0]->nc,
											isochores[0]->ne,
											isochores[0]->UC,
											isochores[0]->DE,
											isochores[0]->md,
											isochores[0]->Md,
											isochores[0]->block);
	  
	  sprintf(mess,"%d Gene Model rules have been read and saved\n",
			  isochores[0]->nclass);
	  printMess(mess);   
	}
  else
	{
	  printMess("Reading Gene Model rules");
	  isochores[0]->nclass = ReadGeneModel(RootFile,
										   isochores[0]->D,
										   isochores[0]->nc,
										   isochores[0]->ne,
										   isochores[0]->UC,
										   isochores[0]->DE,
										   isochores[0]->md,
										   isochores[0]->Md,
										   isochores[0]->block);
	  
	  sprintf(mess,"%d Gene Model rules have been read and saved\n",
			  isochores[0]->nclass);
	  printMess(mess);   
	}

  /* Replication of gene model information for each isochore */
  shareGeneModel(isochores, nIsochores);

  return(nIsochores);
}



