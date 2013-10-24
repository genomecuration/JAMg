/*************************************************************************
*                                                                        *
*   Module: PrintExons                                                   *
*                                                                        *
*   Formatted printing of predicted exons                                *
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

/*  $Id: PrintExons.c,v 1.27 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

extern int GFF;
extern int GFF3;
extern int BP;
extern int PPT;
extern int U12;
extern int SRP;
extern int UTR;
extern float EvidenceFactor;
extern float EvidenceEW;
extern float MRM;



/* Print a predicted exon from a list of exons */
void PrintExon(exonGFF *e, char Name[], char* s, dict* dAA, char* GenePrefix)
{
  char sAux[MAXAA];
  char* rs;
  long p1, p2;
  float kb = 1000.000;
  int nAA = 0;
  char attribute[MAXSTRING] = "";
  char tmpstr[MAXSTRING] = "";
/*   char saux[MAXTYPE]; */
/*   char saux2[MAXTYPE]; */
  
  /* Acquiring real positions of exon */
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;
  
  /* Translation of exon nucleotides into amino acids */
  strcpy(sAux,"\0");
  if((e->Donor->Position >= e->Acceptor->Position)&&(!strcmp(e->Type,sSINGLE)||!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL))){
    if (e->Strand == '+')
      /* Translate codons to amino acids */
      nAA = Translate(p1,p2,
		      e->Frame,
		      (3 - e->Remainder)%3,
		      s, dAA, sAux);
    else
      {
	/* Reverse strand exon */
	if ((rs = (char*) calloc(p2-p1+2,sizeof(char))) == NULL)
	  printError("Not enough memory: translating reverse exon");
      
	/* Reversing and complementing the exon sequence */
	ReverseSubSequence(p1, p2, s, rs);
      
	/* Translate codons to aminoacids */
	nAA = Translate(0,p2-p1,
			e->Frame,
			(3 - e->Remainder)%3,
			rs, dAA, sAux);
	free(rs);
      }
  }
  /* According to selected options formatted output */

  if (GFF3)
    {
      if (! e->evidence){
	if (e->Strand == cFORWARD){
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";start_score=%1.2f;donor_score=%1.2f;donor=%s",e->Acceptor->Score,e->Donor->Score,e->Donor->subtype);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;acceptor_profile_score=%1.2f;acceptor=%s;donor_score=%1.2f;donor=%s",e->Acceptor->Score,e->Acceptor->ScoreAccProfile,e->Acceptor->subtype,e->Donor->Score,e->Donor->subtype);strcat(attribute,tmpstr);
	    if (PPT){
	      sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",e->Acceptor->ScorePPT,e->Acceptor->PositionPPT);strcat(attribute,tmpstr);
	    }
	    if ((BP && e->Donor->class == U2)||(e->Acceptor->class != U2)){
	      sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",e->Acceptor->ScoreBP,e->Acceptor->PositionBP);strcat(attribute,tmpstr);
	    }
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;acceptor_profile_score=%1.2f;acceptor=%s",e->Acceptor->Score,e->Acceptor->ScoreAccProfile,e->Acceptor->subtype);strcat(attribute,tmpstr);
	    if (PPT){
	      sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",e->Acceptor->ScorePPT,e->Acceptor->PositionPPT);strcat(attribute,tmpstr);
	    }
	    if ((BP && e->Donor->class == U2)||(e->Acceptor->class != U2)){
	      sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",e->Acceptor->ScoreBP,e->Acceptor->PositionBP);strcat(attribute,tmpstr);
	    }
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  
	}else{
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";start_score=%1.2f;donor_score=%1.2f;donor=%s",e->Donor->Score,e->Acceptor->Score,e->Acceptor->subtype);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;acceptor_profile_score=%1.2f;acceptor=%s;donor_score=%1.2f;donor=%s",e->Donor->Score,e->Donor->ScoreAccProfile,e->Donor->subtype,e->Acceptor->Score,e->Acceptor->subtype);strcat(attribute,tmpstr);
	    if (PPT){
	      sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",e->Donor->ScorePPT,e->Donor->PositionPPT);strcat(attribute,tmpstr);
	    }
	    if ((BP && e->Donor->class == U2)||(e->Donor->class != U2)){
	      sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",e->Donor->ScoreBP,e->Donor->PositionBP);strcat(attribute,tmpstr);
	    }
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;acceptor_profile_score=%1.2f;acceptor=%s",e->Donor->Score,e->Donor->ScoreAccProfile,e->Donor->subtype);strcat(attribute,tmpstr);
	    if (PPT){
	      sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",e->Donor->ScorePPT,e->Donor->PositionPPT);strcat(attribute,tmpstr);
	    }
	    if ((BP && e->Donor->class == U2)||(e->Donor->class != U2)){
	      sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",e->Donor->ScoreBP,e->Donor->PositionBP);strcat(attribute,tmpstr);
	    }
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  
	}
	sprintf(tmpstr,";coding_potential=%1.2f",e->PartialScore);strcat(attribute,tmpstr);
	if (SRP){sprintf(tmpstr,";homology_score=%1.2f",e->HSPScore);strcat(attribute,tmpstr);}
	if (UTR){sprintf(tmpstr,";rpkm=%1.2f",((e->R)/((p2 - p1 + 1)/kb))/MRM);strcat(attribute,tmpstr);}
      }
      printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=%s_%s_%ld_%ld%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : EXONS,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,
	      e->Type,
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
	      (e->Score==MAXSCORE)? 0.0: e->Score,
	      e->Strand,
	      e->Frame,
	      e->Type,
	      Name,
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
	      attribute);
    }
  else
    {
    if (GFF)
      {
	/* GFF format */
	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s_%s_%ld_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : EXONS,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,
		e->Type,
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0: e->Score,
		e->Strand,
		e->Frame,
		e->Type,
		Name,
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2 );
      }
    else
      {
	/* Default format */
	printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%d\t%s\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		/* saux2,  */
		e->Type,
		(e->evidence)? e->Acceptor->Position: e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		(3 - e->Remainder)%3,
		e->Acceptor->Score,
		e->Donor->Score,
		e->PartialScore,
		e->HSPScore,
		nAA,
		sAux);
      }
  }
}

/* Print a list of exons of the same type */
void PrintExons (exonGFF *e,
                 long ne,
                 int type,
                 char Name[],
                 long l1, long l2,
                 char* Sequence,
                 dict* dAA,
		 char* GenePrefix)
{
  long i;
  char Type[MAXTYPE];
  char strand;
  
  strand = (ne>0)? e->Strand : 'x'; 
  
  Type[0] = '\0';
  switch(type)
    {
    case FIRST: strcpy(Type,sFIRST);
	  break;
    case INTERNAL:strcpy(Type,sINTERNAL);
      break;
    case ZEROLENGTH:strcpy(Type,sZEROLENGTH);
      break;
    case TERMINAL:strcpy(Type,sTERMINAL);
      break;
    case SINGLE:strcpy(Type,sSINGLE);
      break;  
    case ORF:strcpy(Type,sORF);
      break; 
    case UTRFIRST:strcpy(Type,sUTRFIRST);
      break;   
    case UTRFIRSTHALF:strcpy(Type,sUTRFIRSTHALF);
      break;   
    case UTRINTERNAL:strcpy(Type,sUTRINTERNAL);
      break;   
    case UTR5INTERNALHALF:strcpy(Type,sUTR5INTERNALHALF);
      break;   
    case UTR3INTERNALHALF:strcpy(Type,sUTR3INTERNALHALF);
      break;   
    case UTRTERMINALHALF:strcpy(Type,sUTRTERMINALHALF);
      break;   
    case UTRTERMINAL:strcpy(Type,sUTRTERMINAL);
      break;   
    default: strcpy(Type,sEXON);
      strand = 'x'; 
      break;
    }
  
  /* Output header for the following output */
  printf("# %ss(%c) predicted in sequence %s: [%ld,%ld]\n",
		 Type,
		 strand,
		 Name, l1, l2);
  
  /* Print list of exons (type) predicted in the current fragment */
  for (i=0;i<ne;i++)
    PrintExon((e+i), Name, Sequence, dAA,GenePrefix);   
}

/* Print a predicted exon from a assembled gene: gff/geneid format */
void PrintGCDS(exonGFF *e,
                char Name[],
                char* s,
                dict* dAA,
                long ngen,
                int AA1,
                int AA2,
                int nAA,
		int nExon,
		char* GenePrefix
		)
{
/*   float onekb = 1000; */
  long p1,p2;
  float kb=1000.000;
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;

  if (e->Donor->Position == e->Acceptor->Position - 1){return;}
  char attribute[MAXSTRING] = "";
  char tmpstr[MAXSTRING] = "";
  /* Acquiring real positions of exon */
  
/*   char mess[MAXSTRING] = ""; */
/* 	      sprintf(mess,"rpkm: %f",e->R); */
/* 	      printMess(mess); */
  if (GFF3)
    {
      
      /* GFF3 format 5_prime_partial=true ???*/
      if (e->five_prime_partial) { strcpy(attribute,";5_prime_partial=true");}
      if (e->three_prime_partial) { strcat(attribute,";3_prime_partial=true");}
      sprintf(tmpstr,";exon_type=%s",e->Type);strcat(attribute,tmpstr);
      if (! e->evidence){
	if (e->Strand == cFORWARD){
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";start_score=%1.2f;donor_score=%1.2f",e->Acceptor->Score,e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Acceptor->Score,e->Donor->Score);strcat(attribute,tmpstr);
	    
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Acceptor->Score,e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }

	}else{
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";start_score=%1.2f;donor_score=%1.2f",e->Donor->Score,e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Donor->Score,e->Acceptor->Score);strcat(attribute,tmpstr); 
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Donor->Score,e->Acceptor->Score);strcat(attribute,tmpstr); 
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	}
	if (!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)||!strcmp(e->Type,sSINGLE)){
	  sprintf(tmpstr,";coding_potential=%1.2f",e->PartialScore);strcat(attribute,tmpstr);
	}
	if (SRP){sprintf(tmpstr,";homology_score=%1.2f",e->HSPScore);strcat(attribute,tmpstr);}
/* 	if (UTR){sprintf(tmpstr,";rpk=%1.2f",(exp(e->HSPScore/((e->Donor->Position + e->offset2) - (e->Acceptor->Position + e->offset1) +1))-1) * 1000);strcat(attribute,tmpstr);} */
	if (UTR){sprintf(tmpstr,";rpkm=%1.2f",((e->R)/((p2 - p1 + 1)/kb))/MRM);strcat(attribute,tmpstr);}

      }
      /* It's really not necessary to print out exons unless we can predict UTR exons (leave this commented for now) */
/*       printf ("%s\t%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=exon_%s%s_%ld.%i;Parent=mRNA_%s%s_%ld;type=%s\n", */
/* 	      /\* correct stop codon position, Terminal- & Terminal+ *\/  */
/* 	      Name, */
/* 	      (e->evidence)? EVIDENCE : VERSION,      */
/* 	      "exon", */
/* 	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1, */
/* 	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2, */
/* 	      (e->Score==MAXSCORE)? 0.0:e->Score, */
/* 	      e->Strand, */
/* 	      e->Frame, */
/* 	      GenePrefix, */
/* 	      Name, */
/* 	      ngen, */
/* 	      nExon, */
/* 	      GenePrefix, */
/* 	      Name, */
/* 	      ngen, */
/* 	      e->Type); */
      if (!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)||!strcmp(e->Type,sSINGLE)){
      printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=cds_%s%s%s_%ld;Parent=mRNA_%s%s%s_%ld;Target=%s_predicted_protein_%s%s%s_%ld %d %d%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,     
	      "CDS",
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
	      (e->Score==MAXSCORE)? 0.0:e->Score,
	      e->Strand,
	      e->Frame,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      /*nExon,*/
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      VERSION,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      (e->Strand=='+')? nAA-AA2+COFFSET : AA1,
	      (e->Strand=='+')? nAA-AA1+COFFSET : AA2,
	      attribute);
      }else{
	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=utr_%s%s%s_%ld.%i;Parent=mRNA_%s%s%s_%ld%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		"UTR",
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen,
		nExon,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen,
		attribute);
      }
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		e->Type,
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      } else {
	/* Default format for genes */
	printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\tAA %3d:%3d %s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		e->Type,
		(e->evidence)? e->Acceptor->Position :e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		(3 - e->Remainder)%3,
		e->Acceptor->Score,
		e->Donor->Score,
		e->PartialScore,
		e->HSPScore,
		(e->Strand=='+')? nAA-AA2+COFFSET : AA1,
		(e->Strand=='+')? nAA-AA1+COFFSET : AA2,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      }
    }
}
/* Print a predicted exon from a assembled gene: gff/geneid format */
void PrintGUTR(exonGFF *e,
                char Name[],
                char* s,
                long ngen,
		int nExon,
		char* GenePrefix
		)
{
  if (e->Donor->Position == e->Acceptor->Position - 1){return;}
  long p1,p2;
  float kb=1000.000;
  p1 = e->Acceptor->Position + e->offset1 - COFFSET;
  p2 = e->Donor->Position + e->offset2 - COFFSET;

  char attribute[MAXSTRING] = "";
  char tmpstr[MAXSTRING] = "";
/*   char mess[MAXSTRING] = ""; */
  if (GFF3)
    {
      
      /* GFF3 format 5_prime_partial=true ???*/
      if (e->five_prime_partial) { strcpy(attribute,";5_prime_partial=true");}
      if (e->three_prime_partial) { strcat(attribute,";3_prime_partial=true");}
      sprintf(tmpstr,";exon_type=%s",e->Type);strcat(attribute,tmpstr);
      if (! e->evidence){
	if (e->Strand == cFORWARD){
	  
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Acceptor->Score,e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }

	}else{
	  
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",e->Donor->Score,e->Acceptor->Score);strcat(attribute,tmpstr); 
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",e->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",e->Donor->Score);strcat(attribute,tmpstr);
	  }
	}
	if (!strcmp(e->Type,sFIRST)||!strcmp(e->Type,sINTERNAL)||!strcmp(e->Type,sTERMINAL)||!strcmp(e->Type,sSINGLE)){
	  sprintf(tmpstr,";coding_potential=%1.2f",e->PartialScore);strcat(attribute,tmpstr);
	}
	if (SRP){sprintf(tmpstr,";homology_score=%1.2f",e->HSPScore);strcat(attribute,tmpstr);}
	if (UTR){sprintf(tmpstr,";rpkm=%1.2f",((e->R)/((p2 - p1 + 1)/kb))/MRM);strcat(attribute,tmpstr);}

      }

	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=utr_%s%s%s_%ld.%i;Parent=mRNA_%s%s%s_%ld%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		"UTR",
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen,
		nExon,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen,
		attribute);
      
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		e->Type,
		(e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		GenePrefix,
		Name,
		ngen);
      } else {
	/* Default format for genes */
	printf ("%8s %8ld %8ld\t%5.2f\t%c %hd %hd\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		e->Type,
		(e->evidence)? e->Acceptor->Position :e->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		(e->Score==MAXSCORE)? 0.0:e->Score,
		e->Strand,
		e->Frame,
		(3 - e->Remainder)%3,
		e->Acceptor->Score,
		e->Donor->Score,
		e->PartialScore,
		e->HSPScore,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      }
    }
}


/* Print a predicted gene from a assembled gene: gff/geneid format */
void PrintGGene(exonGFF *s,
		exonGFF *e,
                char Name[],
                long ngen,
                float score,
		char* GenePrefix)
{
  if (GFF3)
    {
      /* GFF3 format */
      printf ("%s\t%s%s%s\tgene\t%ld\t%ld\t%.2f\t%c\t.\tID=%s%s%s_%ld;Name=%s%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,     
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? s->Donor->Position : s->Donor->Position + s->offset2,
	      score,
	      e->Strand,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen); 
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s%s%s\tgene\t%ld\t%ld\t%.2f\t%c\t.\t%s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		(e->evidence)? s->Acceptor->Position : s->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		score,
		e->Strand,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      }
    }
}
/* Print a predicted mRNA from a assembled gene: gff/geneid format */
void PrintGmRNA(exonGFF *s,
		exonGFF *e,
                char Name[],
                long ngen,
                float score,
		char* GenePrefix)
{
  if (GFF3)
    {
		
      /* GFF3 format */
      printf ("%s\t%s%s%s\tmRNA\t%ld\t%ld\t%1.2f\t%c\t.\tID=mRNA_%s%s%s_%ld;Parent=%s%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,     
	      (e->evidence)? e->Acceptor->Position : e->Acceptor->Position + e->offset1,
	      (e->evidence)? s->Donor->Position : s->Donor->Position + s->offset2,
	      score,
	      e->Strand,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen);
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s%s%s\tgene\t%ld\t%ld\t%1.2f\t%c\t.\t%s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		(e->evidence)? s->Acceptor->Position : s->Acceptor->Position + e->offset1,
		(e->evidence)? e->Donor->Position : e->Donor->Position + e->offset2,
		score,
		e->Strand,
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      }
    }
}

/* Print a predicted intron from an assembled gene: gff/geneid format */
void PrintGIntron(exonGFF *d,
		  exonGFF *a,
		  char Name[],
		  long ngen,
		  int numInt,
		  char* GenePrefix,
		  int evidence,
		  float score,
		  char* eType)
{
  char attribute[MAXSTRING] = "";
  char tmpstr[MAXSTRING] = "";
  /* float score = 0.0; */
  char intronType[MAXTYPE]; 
  strcpy(intronType,"U2");
  char intronSubtype[MAXTYPE]; 
  strcpy(intronSubtype,"GT-AG");
  short phase = (3 - d->Remainder)%3;
  /* short phase = MIN(0, 3 - a->Frame); */
  long start = (a->evidence)? a->Acceptor->Position: d->Donor->Position + 1 + d->offset2;
  long end = (a->evidence)? a->Donor->Position: a->Acceptor->Position -1 + a->offset1;

  if (evidence){
    /* score = a->Score; */
    
    if(UTR){
      sprintf(tmpstr,";homology_score=%1.2f;rpj=%1.0f",score,((score/EvidenceFactor) - EvidenceEW));strcat(attribute,tmpstr);
      /* sprintf(tmpstr,";homology_score=%1.2f",score);strcat(attribute,tmpstr); */
    }
  }else{
    score = (d->Donor->Score + a->Acceptor->Score)/2;
  }
  if (!(strcmp(d->Donor->subtype,sU12gtag))||!(strcmp(d->Donor->subtype,sU12atac))){
    strcpy(intronType,"U12");
  }
  if (!(strcmp(d->Donor->subtype,sU12gtag))){
    strcpy(intronSubtype,"GT-AG");
  }
  if (!(strcmp(d->Donor->subtype,sU12atac))){
    strcpy(intronSubtype,"AT-AC");
  }
  if (!(strcmp(d->Donor->subtype,sU2gcag))){
    strcpy(intronSubtype,"GC-AG");
  }
  if (d->Strand == '-'){
    phase = (3 - a->Remainder)%3;
  }
  if (GFF3) {
    /* GFF3 format */
      if (! d->evidence && ! a->evidence){
	if(!strcmp(eType,sINTRON)||!strcmp(eType,sUTRINTRON)||!strcmp(eType,sUTR5INTRON)||!strcmp(eType,sUTR3INTRON)){
	  sprintf(tmpstr,";etype=%s",eType);strcat(attribute,tmpstr);
	}
	if (d->Strand == cFORWARD){
	  sprintf(tmpstr,";donor_score=%1.2f;acceptor_score=%1.2f;acceptor_profile_score=%1.2f",d->Donor->Score,a->Acceptor->Score,a->Acceptor->ScoreAccProfile);strcat(attribute,tmpstr);
	    
	  if (PPT){
	    sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",a->Acceptor->ScorePPT,a->Acceptor->PositionPPT);strcat(attribute,tmpstr);
	    }
	  if ((BP && a->Acceptor->class == U2)||(a->Acceptor->class != U2)){
	    sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",a->Acceptor->ScoreBP,a->Acceptor->PositionBP);strcat(attribute,tmpstr);
	  } 
	}else{
	  sprintf(tmpstr,";donor_score=%1.2f;acceptor_score=%1.2f;acceptor_profile_score=%1.2f",a->Acceptor->Score,d->Donor->Score,d->Donor->ScoreAccProfile);strcat(attribute,tmpstr);
	    
	  if (PPT){
	    sprintf(tmpstr,";ppt_score=%1.2f;ppt_pos=%i",d->Donor->ScorePPT,d->Donor->PositionPPT);strcat(attribute,tmpstr);
	    }
	  if ((BP && d->Donor->class == U2)||(d->Donor->class != U2)){
	    sprintf(tmpstr,";bp_score=%1.2f;bp_pos=%i",d->Donor->ScoreBP,d->Donor->PositionBP);strcat(attribute,tmpstr);
	  }
	}
      }
    
    printf ("%s\t%s%s%s\tintron\t%ld\t%ld\t%1.2f\t%c\t%hd\tID=intron_%s%s%s_%ld.%i;Parent=%s%s%s_%ld;type=%s;subtype=%s%s\n",
	    /* correct stop codon position, Terminal- & Terminal+ */ 
	    Name,
	    evidence? EVIDENCE : VERSION,
	    strlen(GenePrefix)>0?"-" : "",
	    GenePrefix,     
	    start,
	    end,
	    (score==MAXSCORE)? 0.0:score,
	    a->Strand,
	    phase,
	    GenePrefix,
	    strlen(GenePrefix)>0?"-" : "",
	    Name,
	    ngen,
	    numInt,
	    GenePrefix,
	    strlen(GenePrefix)>0?"-" : "",
	    Name,
	    ngen,
	    intronType,
	    intronSubtype,
	    attribute
	    );
  } else {
    if (GFF){		 
      /* GFF format */
      printf ("%s\t%s%s%s\t%s_intron\t%ld\t%ld\t%1.2f\t%c\t%hd\t%s%s%s_%ld\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      evidence? EVIDENCE : VERSION,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,     
	      intronType,
	      start,
	      end,
	      (score==MAXSCORE)? 0.0:score,
	      a->Strand,
	      phase,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen);	  	
    }
  }

}
/* Print a predicted intron from an assembled gene: gff/geneid format */
void PrintGExon(  exonGFF *a,
		  int nSegments,
		  char Name[],
		  long ngen,
		  int nExon,
		  char* GenePrefix,
		  int evidence,
		  float score)
{
  exonGFF* d =a;
  exonGFF* e =a;
  char attribute[MAXSTRING] = "";
  char tmpstr[MAXSTRING] = "";
  float kb=1000.000;
  
  if (a->Donor->Position == a->Acceptor->Position - 1){return;}
  float homology_score = 0.0;
  float rpkm = 0.0;
  long start = 0;
  long end = 0;
  float partialScore = 0.0;
  if (nSegments == 1){
    homology_score = (a->HSPScore);
    score = (a->Score);
    partialScore = (a->PartialScore);
    start = (a->evidence)? a->Acceptor->Position: a->Acceptor->Position + a->offset1;
    end = (a->evidence)? a->Donor->Position: a->Donor->Position + a->offset2;
    d = a;
    rpkm = ((a->R)/((end - start + 1)/kb))/MRM; 
  }
  if (nSegments == 2){
    homology_score = (a->PreviousExon->HSPScore + a->HSPScore);
    score = (a->PreviousExon->Score + a->Score);
    partialScore = (a->PreviousExon->PartialScore + a->PartialScore);
    start = (a->PreviousExon->evidence)? a->PreviousExon->Acceptor->Position: a->PreviousExon->Acceptor->Position + a->PreviousExon->offset1;
    end = (a->evidence)? a->Donor->Position: a->Donor->Position + a->offset2;
    d = a->PreviousExon;
    rpkm = ((a->R + a->PreviousExon->R)/((end - start + 1)/kb))/MRM; 
  }
  if (nSegments == 3){
    homology_score = (a->PreviousExon->PreviousExon->HSPScore + a->PreviousExon->HSPScore + a->HSPScore);
    score = (a->PreviousExon->PreviousExon->Score + a->PreviousExon->Score + a->Score);
    partialScore = (a->PreviousExon->PreviousExon->PartialScore + a->PreviousExon->PartialScore + a->PartialScore);
    start = (a->PreviousExon->PreviousExon->evidence)? a->PreviousExon->PreviousExon->Acceptor->Position: a->PreviousExon->PreviousExon->Acceptor->Position + a->PreviousExon->PreviousExon->offset1;
    end = (a->evidence)? a->Donor->Position: a->Donor->Position + a->offset2;
    d = a->PreviousExon->PreviousExon;
    rpkm = ((a->R + a->PreviousExon->PreviousExon->R + a->PreviousExon->R)/((end - start + 1)/kb))/MRM; 
  }
    
  if (GFF3)
    {
      /* GFF3 format 5_prime_partial=true ???*/
      if (e->five_prime_partial) { strcpy(attribute,";5_prime_partial=true");}
      if (e->three_prime_partial) { strcat(attribute,";3_prime_partial=true");}
      if (! e->evidence){
	if (e->Strand == cFORWARD){
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",d->Acceptor->Score,a->Donor->Score);strcat(attribute,tmpstr);
	    
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",d->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",d->Acceptor->Score,a->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Acceptor->Score);strcat(attribute,tmpstr);
	  }

	}else{
	  if (!strcmp(e->Type,sFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",d->Donor->Score,a->Acceptor->Score);strcat(attribute,tmpstr); 
	  }
	  if (!strcmp(e->Type,sTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sSINGLE)){
	    sprintf(tmpstr,";start_score=%1.2f",d->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRFIRST)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRINTERNAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f;donor_score=%1.2f",d->Donor->Score,a->Acceptor->Score);strcat(attribute,tmpstr); 
	  }
	  if (!strcmp(e->Type,sUTR5INTERNALHALF)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Donor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTR3INTERNALHALF)){
	    sprintf(tmpstr,";donor_score=%1.2f",a->Acceptor->Score);strcat(attribute,tmpstr);
	  }
	  if (!strcmp(e->Type,sUTRTERMINAL)){
	    sprintf(tmpstr,";acceptor_score=%1.2f",d->Donor->Score);strcat(attribute,tmpstr);
	  }
	}	
	if (SRP){sprintf(tmpstr,";homology_score=%1.2f",homology_score);strcat(attribute,tmpstr);}
	if(UTR){
	  sprintf(tmpstr,";rpkm=%1.2f",rpkm);strcat(attribute,tmpstr);
	}
      }

      printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%s\tID=exon_%s%s%s_%ld.%i;Parent=mRNA_%s%s%s_%ld%s\n",
	      /* correct stop codon position, Terminal- & Terminal+ */ 
	      Name,
	      (e->evidence)? EVIDENCE : VERSION,
	      strlen(GenePrefix)>0?"-" : "",
	      GenePrefix,     
	      "exon",
	      start,
	      end,
	      (score==MAXSCORE)? 0.0:score,
	      e->Strand,
	      ".",
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      nExon,
	      GenePrefix,
	      strlen(GenePrefix)>0?"-" : "",
	      Name,
	      ngen,
	      attribute);
      
    }
  else
    {
      if (GFF){
	/* GFF format */
	printf ("%s\t%s%s%s\t%s\t%ld\t%ld\t%1.2f\t%c\t%s\t%s%s%s_%ld\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		Name,
		(e->evidence)? EVIDENCE : VERSION,
		strlen(GenePrefix)>0?"-" : "",
		GenePrefix,     
		"exon",
		start,
		end,
		(score==MAXSCORE)? 0.0:score,
		e->Strand,
		".",
		GenePrefix,
		strlen(GenePrefix)>0?"-" : "",
		Name,
		ngen);
      } else {
	/* Default format for genes */
	printf ("%8s %8ld %8ld\t%5.2f\t%c %s %s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n",
		/* correct stop codon position, Terminal- & Terminal+ */ 
		"exon",
		start,
		end,
		(score==MAXSCORE)? 0.0:score,
		e->Strand,
		".",
		".",
		d->Acceptor->Score,
		a->Donor->Score,
		partialScore,
		e->HSPScore);
      }
    }
}

/* Print a predicted exon from a assembled gene: XML format */
void PrintXMLExon(exonGFF *e,
                  char Name[],
                  long ngen,
                  long nExon,
                  int type1,
                  int type2,
		  char* GenePrefix) 
{
  char Type[MAXTYPE];
  
  /* XML format */
  /* exon*/
  printf("      <exon idExon=\"%s%s.G%ldE%ld\" type=\"%s\" frame=\"%hd\" score=\"%.2f\">\n",
		 GenePrefix,Name,ngen,nExon,e->Type,e->Frame,e->Score);
  
  /* Both sites */
  Type[0] = '\0';
  switch(type1)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA: strcpy(Type,sSTA); break;
    case DON: strcpy(Type,sDON); break;
    case STO: strcpy(Type,sSTO); break;
    case POL: strcpy(Type,sPOL); break;
    case TSS: strcpy(Type,sTSS); break;
    case TES: strcpy(Type,sTES); break;
    }
  
  printf("         <site idSite=\"%s%s.G%ldE%ldS1\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 GenePrefix,Name,ngen,nExon,Type,e->Acceptor->Position + e->offset1,e->Acceptor->Score);
  
  Type[0] = '\0';
  switch(type2)
    {
    case ACC: strcpy(Type,sACC); break;
    case STA:strcpy(Type,sSTA); break;
    case DON:strcpy(Type,sDON); break;
    case STO:strcpy(Type,sSTO); break;
    case POL: strcpy(Type,sPOL); break;
    case TSS: strcpy(Type,sTSS); break;
    case TES: strcpy(Type,sTES); break;
    }
  
  printf("         <site idSite=\"%s%s.G%ldE%ldS2\" type=\"%s\" position=\"%ld\" score=\"%.2f\" />\n",
		 GenePrefix,Name,ngen,nExon,Type,e->Donor->Position + e->offset2,e->Donor->Score);
  
  printf("      </exon>\n");
}
