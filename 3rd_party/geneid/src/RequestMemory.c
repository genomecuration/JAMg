/*************************************************************************
*                                                                        *
*   Module: RequestMemory                                                *
*                                                                        *
*   Asking operating system for memory for geneid data structures        *
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

/*  $Id: RequestMemory.c,v 1.19 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Predicted amount of sites and exons found in one split */
extern long NUMSITES, NUMU12SITES, NUMU12EXONS, NUMU12U12EXONS, NUMEXONS, MAXBACKUPSITES, MAXBACKUPEXONS;
extern int scanORF;
extern int SRP,EVD,UTR,GENEID, U12GTAG, U12ATAC;
extern short SPLICECLASSES;

/* Allocating accounting data structure in memory */
account* RequestMemoryAccounting()
{
  account* m; 

  if ((m = (account *) malloc(sizeof(account))) == NULL)
	printError("Not enough memory: account");

  return(m);
}

/* Allocating input sequence in memory */
char* RequestMemorySequence(long L)
{
  char* s;

  if ((s = (char*) calloc(L,sizeof(char))) == NULL)
    printError("Not enough memory: DNA input sequence"); 

  return(s);
}

/* Allocating pack of sites (only in one sense) in memory */
packSites* RequestMemorySites()
{
  packSites* allSites;

  /* Allocating memory for sites */
  if ((allSites = 
       (struct s_packSites *) malloc(sizeof(struct s_packSites))) == NULL)
    printError("Not enough memory: pack of sites");  
  
  /* Start codons */
  if ((allSites->StartCodons = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough memory: start codons");
  
  /* Acceptor sites */
  if ((allSites->AcceptorSites = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough memory: acceptor sites");
	      	      
  /* Donor sites */
  if ((allSites->DonorSites = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough memory: donor sites");

  /* Stop codons */
  if ((allSites->StopCodons = 
       (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
    printError("Not enough memory: stop codons");

  if (UTR){
    /* TSS */
    if ((allSites->TS = 
	 (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
      printError("Not enough memory: TSS");
    
    /* TES */
    if ((allSites->TE = 
	 (struct s_site *) calloc(NUMSITES, sizeof(struct s_site))) == NULL)
      printError("Not enough memory: TES");
  }
  return(allSites);
}

/* Allocating pack of exons (only in one sense) in memory */
packExons* RequestMemoryExons()
{
  packExons* allExons;
  long HowMany;
  
  /* Allocating memory for exons */
  if ((allExons = 
       (struct s_packExons *) malloc(sizeof(struct s_packExons))) == NULL)
    printError("Not enough memory: pack of exons");  

  /* InitialExons */
  HowMany = (long)(NUMEXONS/RFIRST);
  if ((allExons->InitialExons = 
       (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
    printError("Not enough memory: first exons");
      
  /* InternalExons */
  HowMany = (long)(NUMEXONS/RINTER); 
  if ((allExons->InternalExons = 
       (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
    printError("Not enough memory: internal exons");

  /* ZeroLengthExons */
  HowMany = (long)(NUMEXONS/RINTER); 
  if ((allExons->ZeroLengthExons = 
       (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
    printError("Not enough memory: zero length exons");
         
  /* TerminalExons */
  HowMany = (long)(NUMEXONS/RTERMI);
  if ((allExons->TerminalExons = 
       (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
    printError("Not enough memory: terminal exons");

  /* SingleExons */
  HowMany = (long)(NUMEXONS/RSINGL);
  if ((allExons->Singles = 
       (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
    printError("Not enough memory: single genes");
  
  if (UTR){
    /* UTR Exons */
    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->UtrInitialExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: UtrInitialExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->UtrInitialHalfExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: UtrInitialhalfExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->UtrInternalExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: UtrInternalExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->Utr5InternalHalfExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: Utr5InternalHalfExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->Utr3InternalHalfExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: Utr3InternalHalfExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->UtrTerminalHalfExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: UtrTerminalHalfExons");

    HowMany = (long)(NUMEXONS/RUTR);
    if ((allExons->UtrTerminalExons = 
	 (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
      printError("Not enough memory: UtrTerminalExons");
  }
 

  /* IF Scan ORF is switched on... */
  if (scanORF)
    {
	  HowMany = (long)(NUMEXONS/RORF);
	  if ((allExons->ORFs = 
		   (exonGFF*) calloc(HowMany, sizeof(exonGFF))) == NULL)
        printError("Not enough memory: Open Reading Frames");
    }
  allExons->nInitialExons = 0;
  allExons->nInternalExons = 0;
  allExons->nZeroLengthExons = 0;
  allExons->nTerminalExons = 0;
  allExons->nSingles = 0;
  allExons->nORFs = 0;
  allExons->nUtrInitialExons = 0;
  allExons->nUtrInitialHalfExons = 0;
  allExons->nUtrInternalExons = 0;
  allExons->nUtr5InternalHalfExons = 0;
  allExons->nUtr3InternalHalfExons = 0;
  allExons->nUtrTerminalHalfExons = 0;
  allExons->nUtrTerminalExons = 0;

  return(allExons);
}

/* Allocating memory for sorting the set of predicted exons */
exonGFF* RequestMemorySortExons()
{
  exonGFF *exons;
  long HowMany;

  /* Sorting Exons */
  HowMany = NUMEXONS * FSORT;
  if ((exons =
       (exonGFF*) calloc(HowMany, sizeof(exonGFF)))  == NULL)
    printError("Not enough memory: table to sort exons");

  return(exons);
}

/* Allocating memory for sorting the set of predicted sites */
site* RequestMemorySortSites()
{
  site *sites;
  long HowMany;

  /* Sorting Exons */
  HowMany = NUMSITES * FSORT;
  if ((sites =
       (site*) calloc(HowMany, sizeof(site)))  == NULL)
    printError("Not enough memory: table to sort sites");

  return(sites);
}

/* Allocating memory for input evidences (annotations) */
packEvidence* RequestMemoryEvidence()
{
  packEvidence* p;

  /* Allocating memory for structure */
  if ((p = 
       (struct s_packEvidence *) malloc(sizeof(struct s_packEvidence))) == NULL)
    printError("Not enough memory: pack of evidences");    

  /* Evidences sites */
  if ((p->vSites = 
       (struct s_site *) calloc(3*MAXEVIDENCES, sizeof(struct s_site))) == NULL);
//    printError("Not enough memory: evidences sites");

  /* Evidences exons (records) */
  if ((p->vExons =
       (exonGFF*) calloc(MAXEVIDENCES, sizeof(exonGFF)))  == NULL);
//    printError("Not enough memory: evidences exons");

  /* Set counters */
  p->nvSites = 0;
  p->nvExons = 0;

  return(p);
}

HSP* RequestNewHSP()
{
  HSP* p;

  /* New HSP */
  if ((p = (HSP *) malloc(sizeof(HSP))) == NULL)
    printError("Not enough space to hold one new HSP");

  return(p);
}


/* Allocating memory for blast HSPs (homology information) */
packHSP* RequestMemoryHomology()
{
  packHSP* p;
  int i;
  long HowMany;

  /* TWO senses plus THREE reading frames */
  HowMany = STRANDS*FRAMES;

  /* Allocating memory for similarity regions structure */
  if ((p = 
       (struct s_packHSP *) malloc(sizeof(struct s_packHSP))) == NULL)
    printError("Not enough memory: pack of homology information");    

  /* For each (strand,frame) a list of HSPs will be used */
  if ((p->sPairs = 
       (HSP ***) calloc(HowMany, sizeof(HSP**))) == NULL)
    printError("Not enough memory: general array of HSPs");

  for(i=0; i<HowMany; i++)
    if ((p->sPairs[i] = 
		 (HSP **) calloc(MAXHSP, sizeof(HSP*))) == NULL);
// fprintf (stderr, "error with : %i", MAXHSP);
//      printError("Not enough space: individual array of pointers to HSP");  
  
  /* Counters */
  if ((p->nSegments =
       (long*) calloc(HowMany, sizeof(long)))  == NULL)
    printError("Not enough memory: HSP global counter (frame/strand)");

  /* Set counters */
  for(i=0;i<HowMany;i++)
	p->nSegments[i] = 0;
      
  p->nTotalSegments = 0;
  p->visited = 0;

  return(p);
}

/* Alocating memory for external information */
packExternalInformation* RequestMemoryExternalInformation()
{
  packExternalInformation* p;
  int i;
  long HowMany;

  /* TWO senses plus THREE reading frames */
  HowMany = STRANDS*FRAMES;

  /* 0. Allocating main structure */
  /* Allocating memory for similarity regions structure */
  if ((p = 
       (struct s_packExternalInformation *) 
	   malloc(sizeof(struct s_packExternalInformation))) == NULL)
    printError("Not enough memory: pack of external information");  

  /* 1. Dictionary of exon features */
  if ((p->locusNames = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough memory: dictionary of locus names");

  resetDict(p->locusNames);

  /* 2. Arrays of evidence and homology records */
  if (EVD || !GENEID)
	{
	  if ((p->evidence =
		   (packEvidence**) calloc(MAXNSEQUENCES, sizeof(packEvidence*)))  == NULL)
		printError("Not enough memory: array of evidence information");
	  
	  for(i=0; i<MAXNSEQUENCES; i++)
		  p->evidence[i] = (packEvidence*) RequestMemoryEvidence();
	}

  if (SRP)
	{
	  if ((p->homology =
		   (packHSP**) calloc(MAXNSEQUENCES, sizeof(packHSP*)))  == NULL)
		printError("Not enough memory: array of homology information");
	  
	  for(i=0; i<MAXNSEQUENCES; i++)
		p->homology[i] = (packHSP*) RequestMemoryHomology();
	}
  
  /* 3. Counters for every sequence: evidence */
  if (EVD)
	{
	  p->i1vExons = 0;
	  p->i2vExons = 0;
	  p->ivExons = 0;
	}

  /* 4. Counters for every sequence: homology */
  if (SRP)
	{
	  if ((p->iSegments =
		   (long*) calloc(HowMany, sizeof(long)))  == NULL)
		printError("Not enough memory: HSP partial counter (frame/strand)");
	  
	  /* Set counters */
	  for(i=0;i<HowMany;i++)
		p->iSegments[i] = 0;
	  
	  /* Pre-processing array */
	  if ((p->sr = 
		   (float **) calloc(HowMany, sizeof(float*))) == NULL)
		printError("Not enough memory: general preprocessing array of HSPs1");
	  
	  for(i=0; i<HowMany; i++)
		if ((p->sr[i] = 
			 (float *) calloc(LENGTHSi, sizeof(float))) == NULL);
//		  printError("Not enough space: individual preprocessing array of HSPs2"); 
	  if (UTR)
	    {
	      /* Pre-processing array (accurate read counts)*/
	      if ((p->readcount = 
		   (float **) calloc(HowMany, sizeof(float*))) == NULL)
		printError("Not enough memory: general preprocessing array of read counts1");
	  
	      for(i=0; i<HowMany; i++)
		if ((p->readcount[i] = 
		     (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
		  printError("Not enough space: individual preprocessing array of read counts2");
	    }
	}

  return(p);
}


/* Allocating memory for GCinfo */
packGC* RequestMemoryGC()
{
  packGC*  p;

  /* Allocating memory for GC structure */
  if ((p = 
       (struct s_packGC *) malloc(sizeof(struct s_packGC))) == NULL)
    printError("Not enough memory: GC information");    

  /* GC content array */
  if ((p->GC = (long *) calloc(LENGTHSi, sizeof(long))) == NULL)
    printError("Not enough memory: packGC (GC array)");

  /* N's content array */
  if ((p->N = (long *) calloc(LENGTHSi, sizeof(long))) == NULL)
    printError("Not enough memory: packGC (N array)");

  return (p);
}

/* Allocating memory for statistical model parameters of one isochore */
gparam* RequestMemoryParams()
{
  gparam* gp;
  long OligoDim;

  /* 0. Main structure: gparam */
  if ((gp = (gparam*) malloc(sizeof(gparam)))  == NULL)
    printError("Not enough memory: isochore model");

  /* 1. Profiles for signals */
  if ((gp->PolyASignalProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough memory: polyA profile");

  if ((gp->StartProfile = (profile *) malloc(sizeof(profile))) == NULL)  
    printError("Not enough memory: start profile");

  if ((gp->AcceptorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: acceptor profile");

  if ((gp->U12gtagAcceptorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: u12gtag acceptor profile");

  if ((gp->U12atacAcceptorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: u12atac acceptor profile");

  if ((gp->PolyPTractProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: acceptor Poly Pyrimidine Tract profile");

  if ((gp->BranchPointProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: acceptor Branch Point profile");

  if ((gp->U12BranchPointProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: u12 acceptor Branch Point profile");
	
  if ((gp->DonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 donor profile");

  if ((gp->U2gcagDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 GCAG donor profile");
	
  if ((gp->U2gcagDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 GCAG donor profile");
	
  if ((gp->U2gtaDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 GTA donor profile");
	
  if ((gp->U2gtgDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 GTG donor profile");
	
  if ((gp->U2gtyDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U2 GTY donor profile");

  if ((gp->U12gtagDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U12 GTAG donor profile");

  if ((gp->U12atacDonorProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: U12 ATAC donor profile");

  if ((gp->StopProfile = (profile *) malloc(sizeof(profile))) == NULL)
    printError("Not enough memory: stop profile");

  /* 2. Markov model: initial and transition values */
  OligoDim = (int)pow((float)4,(float)OLIGOLENGTH);
  
  if ((gp->OligoLogsIni[0]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov initial pentanucleotides (0)");

  if ((gp->OligoLogsIni[1]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov initial pentanucleotides (1)");

  if ((gp->OligoLogsIni[2]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov initial pentanucleotides (2)");

  if ((gp->OligoLogsTran[0]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov chains - hexanucleotides (0)");

  if ((gp->OligoLogsTran[1]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov chains - hexanucleotides (1)");

  if ((gp->OligoLogsTran[2]=(float *) calloc(OligoDim, sizeof(float))) == NULL)
    printError("Not enough memory: Markov chains - hexanucleotides (2)"); 

  /* 3. Markov temporary data structures to compute every split: LENGTHSi */
  if ((gp->OligoDistIni[0] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. initial pentanucleotides sum (0)");
  
  if ((gp->OligoDistIni[1] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. initial pentanucleotides sum (1)");
  
  if ((gp->OligoDistIni[2] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. initial pentanucleotides sum (2)");
  
  if ((gp->OligoDistTran[0] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. Markov hexanucleotides sum (0)");    
  
  if ((gp->OligoDistTran[1] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. Markov hexanucleotides sum (1)");    
  
  if ((gp->OligoDistTran[2] = 
       (float *) calloc(LENGTHSi, sizeof(float))) == NULL)
    printError("Not enough memory: temp. Markov hexanucleotides sum (2)");    

  /* 4. Exons score parameters */
  if ((gp->Initial = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough memory: exons scoring parameters (first)");

  if ((gp->Internal = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough memory: exons scoring parameters (internal)");
 
  if ((gp->Terminal = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough memory: exons scoring parameters (terminal)");

  if ((gp->Single = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough memory: exons scoring parameters (single)");
  
  if ((gp->utr = (paramexons *)malloc(sizeof(paramexons))) == NULL)  
    printError("Not enough memory: exons scoring parameters (utr)");

  return(gp);
}

/* Allocating memory (all of the isochores) */
gparam ** RequestMemoryIsochoresParams()
{
  gparam** isochores;
  int i;

  /* Allocating the array of isochores */
  if ((isochores = (gparam **) calloc(MAXISOCHORES, sizeof(gparam *))) == NULL)
    printError("Not enough memory: isochores array");

  /* Allocating every separate isochore */
  for(i=0; i<MAXISOCHORES; i++)
    isochores[i] = (gparam *) RequestMemoryParams();
  
  /* Allocating space for global parameters (gene model) */
  if ((isochores[0]->D = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough memory: dictionary of exon types");

  if ((isochores[0]->nc = (int *) calloc(MAXENTRY, sizeof(int))) == NULL)
    printError("Not enough memory: nc-array (gene model)");

  if ((isochores[0]->ne = (int *) calloc(MAXENTRY, sizeof(int))) == NULL)
    printError("Not enough memory: ne-array (gene model)");

  if ((isochores[0]->md = (long *) calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough memory: minDist (gene model)");

  if ((isochores[0]->Md = (long *) calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough memory: maxDist (gene model)");

  return(isochores);
}

/* Allocating memory for Markov chain in every position of the profile */
/* This function is called from readparams.c when dimension is known */
void RequestMemoryProfile(profile* p)
{
  int i;
  
  /* Transition probabilities in every position of the PWA */
  p->dimensionTrans = (int)pow((float)5,(float)(p->order+1));   

  for(i=0; i < p->dimension; i++)
    if ((p->transitionValues[i] = 
		 (float *) calloc(p->dimensionTrans, sizeof(float))) == NULL)
      printError("Not enough memory: signal profile transition values");
}

/* Allocating memory for the best set of predicted genes and extra info */
packGenes* RequestMemoryGenes()
{
  packGenes* pg;
  int aux, aux2, aux3;

  /* 0. Allocating memory for pack of genes (main structure) */
  if ((pg = 
       (struct s_packGenes *) malloc(sizeof(struct s_packGenes))) == NULL)
    printError("Not enough memory: pack of genes");  

  /* 1. Allocating memory space for Ghost Exon */
  if (( pg->Ghost = (exonGFF *) malloc(sizeof(exonGFF))) == NULL)
    printError("Not enough memory: Ghost Exon");

  if ((pg->Ghost->Acceptor = (site *) malloc(sizeof(site))) == NULL)
    printError("Not enough memory: Ghost Exon acceptor");

  if ((pg->Ghost->Donor = (site *) malloc(sizeof(site))) == NULL)
    printError("Not enough memory: Ghost Exon donor");

  /* Mark this exon as Ghost Exon */
  pg->Ghost->Strand = '*';
  pg->Ghost->Remainder = 0;
  pg->Ghost->Frame = 0;
  pg->Ghost->GeneScore = 0.0;
  pg->Ghost->Donor->Score = 0;
  pg->Ghost->Donor->Position = 0;
  pg->Ghost->Donor->class = U2;
  pg->Ghost->Acceptor->Score = 0;
  pg->Ghost->Acceptor->Position = 0;
  pg->Ghost->Acceptor->class = U2;
  strcpy(pg->Ghost->Donor->type,sU2type);
  strcpy(pg->Ghost->Donor->subtype,sU2);
  strcpy(pg->Ghost->Acceptor->type,sU2type);
  strcpy(pg->Ghost->Acceptor->subtype,sU2);
  pg->Ghost->offset1 = 0;
  pg->Ghost->offset2 = 0;
  pg->GOptim = pg->Ghost;
  
  /* 2. Allocating memory space for Ga */
  /* Ga is the array of best predicted genes (in every gene class) */ 
  if ((pg->Ga = (exonGFF* ***)calloc(MAXENTRY, sizeof(exonGFF* **))) == NULL)
    printError("Not enough memory: Ga array of genes");

  /* Initialize Ga-exons: everybody looking at the Ghost exon */
  /* MAXENTRY represents the maximum number of gene classes */
  for(aux=0; aux<MAXENTRY; aux++)
    {
      if ((pg->Ga[aux] = (exonGFF* **)calloc(FRAMES, sizeof(exonGFF* *))) == NULL)
        printError("Not enough memory: 6 frames in Ga array of genes");

      for(aux2=0; aux2 < FRAMES; aux2++){
      if ((pg->Ga[aux][aux2] = (exonGFF* *)calloc(SPLICECLASSES, sizeof(exonGFF*))) == NULL)
        printError("Not enough memory: 3 splice classes in Ga array of genes");

	for(aux3=0; aux3 < SPLICECLASSES; aux3++){
	  pg->Ga[aux][aux2][aux3] = pg->Ghost;
	}
      }
      	
    }

  /* 3. Allocate memory space for the set of auxiliary arrays */
  /* Memory for the array of sorting by donor functions (one per class) */
  if ((pg->d = (exonGFF* **)calloc(MAXENTRY, sizeof(exonGFF* *))) == NULL)
    printError("Not enough memory: set of d-arrays (sort by donor)");
  
  /* Memory for every sorting function (alone) */
  for(aux=0; aux < MAXENTRY; aux++) 
    if ((pg->d[aux] = (exonGFF* *)calloc(FDARRAY * NUMEXONS, sizeof(exonGFF*))) == NULL)
      printError("Not enough memory: sort-by-donor functions");
  
  if ((pg->km = (long *)calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough memory: total counters of sort-by-donor functions");

  if ((pg->je = (long *)calloc(MAXENTRY, sizeof(long))) == NULL)
    printError("Not enough memory: partial counters of sort-by-donor functions");

  return(pg);
}

/* Allocating memory for temporary but necessary backup information */
/* This information is required to continue the process between 2 splits */
packDump* RequestMemoryDumpster()
{
  packDump* d;
  long HowMany;

  /* 0. Allocating memory for dumpster */
  if ((d = 
       (struct s_packDump *) malloc(sizeof(struct s_packDump))) == NULL)
    printError("Not enough memory: dumpster");  

  /* 1. Temporary dumpster Sites */
  if ((d->dumpSites = 
       (struct s_site *) calloc(MAXBACKUPSITES, sizeof(struct s_site))) == NULL)
;//    printError("Not enough memory: backup sites");

  /* 2. Temporary dumpster exons */
  if ((d->dumpExons = 
       (exonGFF*) calloc(MAXBACKUPEXONS, sizeof(struct s_exonGFF))) == NULL)
;//    printError("Not enough memory: backup exons");  

  /* 3. Dumpster hash to find backup exons quickly */
  if ((d->h = 
       (dumpHash*) malloc(sizeof(struct s_dumpHash))) == NULL)
    printError("Not enough memory: dumpster hash table structure");  
  
  HowMany = (long) (MAXBACKUPEXONS / HASHFACTOR);

  if ((d->h->T =
       (dumpNode**) calloc(HowMany, sizeof(dumpNode*))) == NULL)
;//    printError("Not enough memory: dumpster hash table");   

  /* 4. Set counters */
  d->ndumpSites = 0;
  d->ndumpExons = 0;
  resetDumpHash(d->h);
  
  return(d);
}

/* Allocate the genetic code to translate from genes into proteins */
dict* RequestMemoryAaDictionary()
{
  dict* dAA;

  /* 1. Allocating memory for aminoacid dictionary */
  if ((dAA = (dict *)malloc(sizeof(dict))) == NULL)
    printError("Not enough memory: aa-dictionary");

  /* 2. Reset dictionary */
  resetDict(dAA);

  /* 3. Filling it in: amino acids and codons */
  setAADict(dAA,"GCA",'A');
  setAADict(dAA,"GCC",'A');
  setAADict(dAA,"GCG",'A');
  setAADict(dAA,"GCT",'A');
  setAADict(dAA,"TGC",'C');
  setAADict(dAA,"TGT",'C');
  setAADict(dAA,"GAC",'D');
  setAADict(dAA,"GAT",'D');
  setAADict(dAA,"GAA",'E');
  setAADict(dAA,"GAG",'E');
  setAADict(dAA,"TTC",'F');
  setAADict(dAA,"TTT",'F');
  setAADict(dAA,"GGA",'G');
  setAADict(dAA,"GGC",'G');
  setAADict(dAA,"GGG",'G');
  setAADict(dAA,"GGT",'G');
  setAADict(dAA,"CAC",'H');
  setAADict(dAA,"CAT",'H');
  setAADict(dAA,"ATA",'I');
  setAADict(dAA,"ATC",'I');
  setAADict(dAA,"ATT",'I');
  setAADict(dAA,"AAA",'K');
  setAADict(dAA,"AAG",'K');
  setAADict(dAA,"TTA",'L');
  setAADict(dAA,"TTG",'L');
  setAADict(dAA,"CTA",'L');
  setAADict(dAA,"CTC",'L');
  setAADict(dAA,"CTG",'L');
  setAADict(dAA,"CTT",'L');
  setAADict(dAA,"ATG",'M');
  setAADict(dAA,"AAC",'N');
  setAADict(dAA,"AAT",'N');
  setAADict(dAA,"CCA",'P');
  setAADict(dAA,"CCC",'P');
  setAADict(dAA,"CCG",'P');
  setAADict(dAA,"CCT",'P');
  setAADict(dAA,"CAA",'Q');
  setAADict(dAA,"CAG",'Q');
  setAADict(dAA,"AGA",'R');
  setAADict(dAA,"AGG",'R');
  setAADict(dAA,"CGA",'R');
  setAADict(dAA,"CGC",'R');
  setAADict(dAA,"CGG",'R');
  setAADict(dAA,"CGT",'R');
  setAADict(dAA,"AGC",'S');
  setAADict(dAA,"AGT",'S');
  setAADict(dAA,"TCA",'S');
  setAADict(dAA,"TCC",'S');
  setAADict(dAA,"TCG",'S');
  setAADict(dAA,"TCT",'S');
  setAADict(dAA,"ACA",'T');
  setAADict(dAA,"ACC",'T');
  setAADict(dAA,"ACG",'T');
  setAADict(dAA,"ACT",'T');
  setAADict(dAA,"GTA",'V');
  setAADict(dAA,"GTC",'V');
  setAADict(dAA,"GTG",'V');
  setAADict(dAA,"GTT",'V');
  setAADict(dAA,"TGG",'W');
  setAADict(dAA,"TAC",'Y');
  setAADict(dAA,"TAT",'Y');
  setAADict(dAA,"TAA",'*');
  setAADict(dAA,"TAG",'*');
  setAADict(dAA,"TGA",'*');

  return(dAA);
}


