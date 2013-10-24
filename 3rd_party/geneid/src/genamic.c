/*************************************************************************
*                                                                        *
*   Module: genamic                                                      *
*                                                                        *
*   Assembling genes from the input set of exons                         *
*                                                                        *
*   This file is part of the geneid 1.4 Distribution                     *
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

/*  $Id: genamic.c,v 1.20 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Complete gene prediction (sites and exons) or only assembling */
extern int GENEID;
extern int RSS;
extern float U12_SPLICE_SCORE_THRESH;
extern float U12_EXON_SCORE_THRESH;

void genamic(exonGFF* E, long nExons, packGenes* pg, gparam* gp)
{
  long i,j,j2;
  short frame,remainder,spliceclass,dclass;
  int h;
  char saux[MAXTYPE];
  int type,etype;
  long MaxDist;
  long MinDist;
  char mess[MAXSTRING];
  int current_exon_is_u12 = 0;
  int thresholdmet = 1;

  /* 0. Starting process ... */
  printMess("-- Running gene assembling (genamic) --");

  /* geneid sends this set of exonsGFF together with the Ga and d-array*/
  sprintf(mess,"%ld exons GFF received from geneid",nExons);
  printMess(mess);

  /* Frame/remainder of reverse-strand exons must be exchanged */
  if (GENEID)
    {
      /* Exons from the current fragment of sequence */
      SwitchFrames(E,nExons);
	  
      /* Exons from the last fragment of sequence */
      SwitchFramesDa(pg,gp->nclass);
    }
  else
    {
      /* GENAMIC processing only: option -O */
      /* Exons from the current fragment of sequence */
      SwitchFrames(E,nExons);
    }
  
  /* 1. Create a set of sorted by donor array of pointer to exons */
  printMess("Sorting exons by donor...");
  
  /* Build a set of sorting exons by donor functions */
  BuildSort(gp->D, gp->nc, gp->ne, 
	    gp->UC, gp->DE, gp->nclass, 
	    pg->km, pg->d, E, nExons);
  
  /* 2. Genamic Algorithm in linear time (size of input) */
  printMess("Assembling Genes...");

  /* For every exon verify non-blocking and distances well defined */
  /* and after this, it might be assemble with the best gene computed */
  for(i=0; i< nExons; i++)
    {
      /* What is the type of this exon? */
      if (strcmp((E+i)->Type,"")){
	saux[0]='\0';
	strcpy (saux, (E+i)->Type);
	strcat (saux, &((E+i)->Strand));
	type = getkeyDict(gp->D,saux);
	frame = (E+i)->Frame;      
	(E+i)->PreviousExon = pg->Ghost;
	(E+i)->GeneScore = (E+i)->Score;
	current_exon_is_u12 = 0;
	spliceclass = (E+i)->Acceptor->class;
	if (!strcmp((E+i)->Type,sEND) || !strcmp((E+i)->Type,sBEGIN)|| !strcmp((E+i)->Type,sSINGLE)){
	  current_exon_is_u12 = 0;
	}else{ 
	  if ((E+i)->Acceptor->class == U2){
	    current_exon_is_u12 = 0;
	  } else { 
	    if (((E+i)->Acceptor->class == U12gtag)||((E+i)->Acceptor->class == U12atac)){
	      current_exon_is_u12 = 1;
	    }
	  }
	}
	if (type != NOTFOUND)
	  {
	    /* For every equivalent class building the best gene ending with it */
	    for(h=0; h < gp->ne[type]; h++)
	      {
		etype = gp->DE[type][h];
		j = pg->je[etype];
		MaxDist = gp->Md[etype];
		MinDist = gp->md[etype];
		thresholdmet = 1;			  

		/* Checking maximum distance allowed requirement */
		if ((MaxDist != INFI) &&
		    (pg->Ga[etype][frame][spliceclass]->Strand !='*') &&
		    (pg->Ga[etype][frame][spliceclass]->Donor->Position 
		     + 
		     pg->Ga[etype][frame][spliceclass]->offset2) 
		    < 
		    ((E+i)->Acceptor->Position + (E+i)->offset1)
		    + ((E+i)->evidence - pg->Ga[etype][frame][spliceclass]->evidence)
		    - MaxDist) 
		  {
		    /* loop backward searching another best gene matching MAX distance */
		    pg->Ga[etype][frame][spliceclass] = pg->Ghost; 
		    j2 = j-1;
		    while (j2>=0 && j2 < pg->km[etype]  && 
			   ((pg->d[etype][j2]->Donor->Position 
			     + pg->d[etype][j2]->offset2)
			    >= 
			    ((E+i)->Acceptor->Position + (E+i)->offset1)
			    + ((E+i)->evidence - pg->d[etype][j2]->evidence)
			    - MaxDist))
		      {
			remainder = pg->d[etype][j2]->Remainder;
			dclass = pg->d[etype][j2]->Donor->class;

			if ((pg->d[etype][j2]->GeneScore > 
			     pg->Ga[etype][remainder][dclass] -> GeneScore)
			    ){
			  pg->Ga[etype][remainder][dclass] = pg->d[etype][j2];
			}
			j2--;
		      }
		  }
			  
		/* Loop forward: One scan over each donor-sort array */
		/* while minimum distance allowed requirement is OK */
		/* Update best partial genes between current and previous exon */
		while(j < pg->km[etype] &&
		      ((pg->d[etype][j]->Donor->Position 
			+ pg->d[etype][j]->offset2)
		       <= 
		       ((E+i)->Acceptor->Position + (E+i)->offset1) 
		       + ((E+i)->evidence - pg->d[etype][j]->evidence)
		       - MinDist)
		      )
		  {
		    remainder = pg->d[etype][j]->Remainder;
		    dclass = pg->d[etype][j]->Donor->class;
		    if ((frame == remainder && spliceclass == dclass &&
			 ((pg->d[etype][j]->Donor->Position 
			   + pg->d[etype][j]->offset2)
			  < 
			  ((E+i)->Acceptor->Position + (E+i)->offset1)
			  + ((E+i)->evidence - pg->d[etype][j]->evidence) 
			  - MaxDist))
			)
		      {
			/* Skip this exon because max distance not ok */
		      }
		    else
		      {
			if (pg->d[etype][j]->GeneScore > pg->Ga[etype][remainder][dclass]->GeneScore) 
			  {
			    pg->Ga[etype][remainder][dclass] = pg->d[etype][j];
			  }
		      }
		    j++;
		  }
		pg->je[etype] = j;

		/* Assembling the exon with the best compatible gene before it */
		/* Verify group rules if there are evidence exons (annotations) */
		if (current_exon_is_u12){
		  if (pg->Ga[etype][frame][spliceclass]->Donor->class != U2){				  
		    if((((pg->Ga[etype][frame][spliceclass]->Donor->Score + (E+i)->Acceptor->Score) > U12_SPLICE_SCORE_THRESH)
			&&
			((pg->Ga[etype][frame][spliceclass]->Score + (E+i)->Score) > U12_EXON_SCORE_THRESH)
			)
		       ||
		       (
			(E+i)->evidence || pg->Ga[etype][frame][spliceclass]->evidence
			)
		       ){
		      thresholdmet = 1;
		    } else {
		      thresholdmet = 0;
		    }
		  } else {			  
		    thresholdmet = 0;
		  }			  	
		} else {
		  if (pg->Ga[etype][frame][spliceclass]->Donor->class != U2){
		    thresholdmet = 0;
		  }
		}

		if ((!(strcmp(pg->Ga[etype][frame][spliceclass]->Group,(E+i)->Group))
		     || gp->block[etype] == NONBLOCK)
		    &&
		    ((pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score) > (E+i)->GeneScore)
		    &&
		    (thresholdmet) 
		    )
		  {
			  
		    if (!strcmp((E+i)->Type,sINTRON)||!strcmp((E+i)->Type,sUTRINTRON)||!strcmp((E+i)->Type,sUTR5INTRON)||!strcmp((E+i)->Type,sUTR3INTRON)){			    
		      (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
		      (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
		      (E+i)->lValue = pg->Ga[etype][frame][spliceclass]->lValue;
		      (E+i)->rValue = pg->Ga[etype][frame][spliceclass]->rValue;
		    }else
		      {if (RSS && ((E+i)->Donor->Position == (E+i)->Acceptor->Position -1)){
			  (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
			  (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
			  (E+i)->Frame = pg->Ga[etype][frame][spliceclass]->Frame;
			  (E+i)->Remainder = pg->Ga[etype][frame][spliceclass]->Remainder;
			  (E+i)->lValue = pg->Ga[etype][frame][spliceclass]->lValue;
			  (E+i)->rValue = pg->Ga[etype][frame][spliceclass]->rValue;
			}else{

			 if (((E+i)->Strand == '+') && 
			     ((pg->Ga[etype][frame][spliceclass]->rValue == 1 && (E+i)->lValue == 1)
			      ||
			      (pg->Ga[etype][frame][spliceclass]->rValue == 2 && (E+i)->lValue == 2)
			      ||
			      (pg->Ga[etype][frame][spliceclass]->rValue == 3 && ((E+i)->lValue == 2 || (E+i)->lValue == 3))))
			   {
			     /* FWD: Avoiding building a stop codon */  
			   }
			 else
			   {
			     if (((E+i)->Strand == '-') && 
				 ((pg->Ga[etype][frame][spliceclass]->lValue == 1 && (E+i)->rValue == 1)
				  ||
				  (pg->Ga[etype][frame][spliceclass]->lValue == 2 && (E+i)->rValue == 2)
				  ||
				  ((pg->Ga[etype][frame][spliceclass]->lValue == 2 || pg->Ga[etype][frame][spliceclass]->lValue == 3) && (E+i)->rValue == 3)))
			       {
				 /* RVS: Avoiding building a stop codon */
			       }
			     else
			       {
				 (E+i)->GeneScore = pg->Ga[etype][frame][spliceclass]->GeneScore + (E+i)->Score;
				 (E+i)->PreviousExon = pg->Ga[etype][frame][spliceclass];
			       }
			   }
		       }
		      }
		  }
	      }

	    /* Updating the best gene assembled (final gene) */
	    if ((((E+i)->GeneScore) > (pg->GOptim -> GeneScore))){
	      if (((E+i)->PreviousExon->Strand == '*')&&(!strcmp((E+i)->Type,sINTRON) || !strcmp((E+i)->Type,sUTRINTRON) || !strcmp((E+i)->Type,sUTR5INTRON) || !strcmp((E+i)->Type,sUTR3INTRON))){
	      }else{
		pg->GOptim = (E+i);
	      }
	    }
	  }
      }
    }

  /* 3. Undo the change between frame and remainder */
  if (GENEID)
    {
      /* Exons predicted in the current fragment of sequence */
      SwitchFrames(E,nExons);
	  
      /* Only changing exons predicted in the last fragment of sequence */
      SwitchFramesDb(pg,gp->nclass);
    }
  else
    {
      /* GENAMIC processing only: option -O */
      UndoFrames(E,nExons);
    }
  
  /* Finishing process */
  printMess("-- Finishing gene assembling (genamic) --\n");
}

