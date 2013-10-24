/*************************************************************************
*                                                                        *
*   Module: beggar                                                       *
*                                                                        *
*   Estimate the memory needed to run geneid (current configuration)     *
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

/*  $Id: beggar.c,v 1.7 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

extern long NUMSITES, NUMU12SITES, NUMEXONS, NUMU12EXONS, NUMU12U12EXONS, MAXBACKUPSITES, MAXBACKUPEXONS;
extern int U12;
/* Computing the memory required to execute geneid */
void beggar(long L)
{
  float memTotal;
  float memSites;
  float memExons;
  float memEvi;
  float memHSP;
  float memGC;
  float memGenes;
  float memParams;
  float memSequence;
  float memBackup;
  int numprofiles = 16;
  int spliceclasses = 1;
  if (U12){
    spliceclasses = 3;
  }
  /* packSites: +/- */
  memSites = STRANDS * (sizeof(struct s_packSites) + (4 * NUMSITES * sizeof(struct s_site)));
  
  /* packExons: +/- */
  memExons = STRANDS * (sizeof(struct s_packExons) +
						((NUMEXONS/RFIRST + NUMEXONS/RINTER +
						  NUMEXONS/RTERMI + NUMEXONS/RSINGL + NUMEXONS/RUTR) *
						 sizeof(exonGFF))
			);
  /* Sort exons */
  memExons += NUMEXONS * FSORT * sizeof(exonGFF);
  
  /* Evidences */
  memEvi = (MAXEVIDENCES * sizeof(exonGFF)) +
	(3*MAXEVIDENCES * sizeof(struct s_site));
  
  /* HSPs */
  memHSP =  sizeof(struct s_packHSP) +
	(STRANDS * FRAMES * MAXHSP * sizeof(HSP));
  
  /* G+C INFO */
  memGC = 2 * LENGTHSi * sizeof(long);
  
  /* pack Genes: ghost, Ga, d, km-jm */
  memGenes = sizeof(struct s_packGenes) +
	2 * (sizeof(exonGFF) + 2 * sizeof(struct s_site)) +
	MAXENTRY * FRAMES * spliceclasses * sizeof(exonGFF*) +
	MAXENTRY * (2 * NUMEXONS) * sizeof(exonGFF*) +
	2 * MAXENTRY * sizeof(long);
  
  /* Statistical model: profiles, Markov model, Markov tmp, exon values,  */
  /* isochores and gene model */
  memParams = (numprofiles * (sizeof(profile) + AVG_DIM * AVG_ORDER * sizeof(float))) 
	+ (6 * OLIGO_DIM * sizeof(float)) 
	+ (6 * LENGTHSi * sizeof(float)) 
	+ (4 * sizeof(paramexons));
  
  memParams += MAXISOCHORES * sizeof(gparam *);
  memParams = memParams * MAXISOCHORES;

  memParams += sizeof(dict) + 2*MAXENTRY* sizeof(int) +
	2*MAXENTRY* sizeof(long);
            
  /* Backup operations */
  memBackup = sizeof(struct s_packDump) + MAXBACKUPSITES * sizeof(struct s_site)
	+ MAXBACKUPEXONS * sizeof(exonGFF) + sizeof(struct s_dumpHash)
	+ ((MAXBACKUPEXONS / HASHFACTOR) * sizeof(dumpNode*));

  /* Sequence space */
  memSequence = L * sizeof(char);


  memTotal = 
	memSites +
	memExons +
	memEvi +
	memHSP +
	memGC +
	memGenes +
	memParams +
	memSequence +
	memBackup;

  
  /* Display numbers */
  printf("AMOUNT of MEMORY required by current geneid configuration\n");
  printf("---------------------------------------------------------\n\n");
  
  printf("Sites\t\t\t: %.2f Mb\n",(float)memSites/(float)MEGABYTE);
  printf("Exons\t\t\t: %.2f Mb\n",(float)memExons/(float)MEGABYTE);
  printf("Evidences\t\t: %.2f Mb\n",(float)memEvi/(float)MEGABYTE);
  printf("Homology\t\t: %.2f Mb\n",(float)memHSP/(float)MEGABYTE);
  printf("G+C info\t\t: %.2f Mb\n",(float)memGC/(float)MEGABYTE);
  printf("Genes\t\t\t: %.2f Mb\n",(float)memGenes/(float)MEGABYTE);
  printf("Statistical model\t: %.2f Mb\n",(float)memParams/(float)MEGABYTE);
  printf("Backup operations\t: %.2f Mb\n",(float)memBackup/(float)MEGABYTE);
  printf("Input sequence\t\t: %.2f Mb\n",(float)memSequence/(float)MEGABYTE);

  printf("---------------------------------------------------------\n\n");
  printf("TOTAL AMOUNT\t\t: %.2f Mb\n",
		 memTotal/MEGABYTE);

  exit(0);
}
