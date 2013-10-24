/*************************************************************************
*                                                                        *
*   Module: SwitchFrames                                                 *
*                                                                        *
*   Exchange frame and remainder from reverse exons                      *
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

/*  $Id: SwitchFrames.c,v 1.9 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Exchange frame and remainder in the input exons */
void SwitchFrames(exonGFF* e, long n)
{
  long i;
  int f;
/*   short class; */
/*   float score; */
 
  /* Exchange frame/rmd in reverse exons and reset the selected flag */
  for (i=0; i<n; i++)
    {      
      if ((e+i)->Strand == '-')
	{
	  if (!strcmp((e+i)->Type,sINTRON)||!strcmp((e+i)->Type,sUTR5INTRON)||!strcmp((e+i)->Type,sUTR3INTRON)||!strcmp((e+i)->Type,sUTRINTRON)){
	    f=(e+i)->Frame;
	    (e+i)->Remainder=f;
	    /* (e+i)->Frame = ( 3 - (e+i)->Remainder )%3; */
	  }else{
	    f=(e+i)->Frame;
	    (e+i)->Frame=(e+i)->Remainder;
	    (e+i)->Remainder=f;
	  }
/* 	  score=(e+i)->Donor->Score; */
/* 	  class=(e+i)->Donor->class; */
/* 	  (e+i)->Donor->Score = (e+i)->Acceptor->Score; */
/* 	  (e+i)->Donor->class = (e+i)->Acceptor->class; */
/* 	  (e+i)->Acceptor->Score = score; */
/* 	  (e+i)->Acceptor->class = class; */
	}else{
	if (!strcmp((e+i)->Type,sINTRON)||!strcmp((e+i)->Type,sUTR5INTRON)||!strcmp((e+i)->Type,sUTR3INTRON)||!strcmp((e+i)->Type,sUTRINTRON)){
	    f=(e+i)->Frame;
	    (e+i)->Remainder=f;
	  }
      }

      /* Mark exon as prediction in the current fragment */
      (e+i)->selected = 0;
    }

}

/* Exchange frame and remainder in the sorted-by-donor exons only once */
/* Right now, d-array only contain exons from last fragment processing */
void SwitchFramesDa(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;
/*   short class; */
/*   float score; */

  /* Screening every class looking for exons... */
  for (i=0; i < nclass; i++) 
    {
      /* Traversing the list of exons in this class */
      for (j=0; j < pg->km[i]; j++)
	if (pg->d[i][j]->Strand == '-')
	  {
            /* Exchange frame/rmd only once */
            /* One exon might be in more than one list */
            if (!pg->d[i][j]->selected)
	      {	    
		f= pg->d[i][j]->Frame;
		pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
		pg->d[i][j]->Remainder = f;
/* 		score=pg->d[i][j]->Donor->Score; */
/* 		class=pg->d[i][j]->Donor->class; */
/* 		pg->d[i][j]->Donor->Score = pg->d[i][j]->Acceptor->Score; */
/* 		pg->d[i][j]->Donor->class = pg->d[i][j]->Acceptor->class; */
/* 		pg->d[i][j]->Acceptor->Score = score; */
/* 		pg->d[i][j]->Acceptor->class = class; */
		/* Mark exon */
		pg->d[i][j]->selected = 1;
	      }
	  }
    }
 
}   

/* Restore original frame and remainder in exons from last fragment */
void SwitchFramesDb(packGenes* pg, int nclass)
{
  long i;
  long j;
  int f;
/*   short class; */
/*   float score; */

  /* Screening every class looking for exons... */
  for (i=0; i < nclass; i++)
    {
      /* Traversing the list of exons in this class */
      for (j=0; j < pg->km[i]; j++)
	if (pg->d[i][j]->Strand == '-')
	  {
	    /* Exchange frame/rmd only once */
	    /* Only exons from last fragment will have selected = 1 */
	    if (pg->d[i][j]->selected)
	      {
		f = pg->d[i][j]->Frame;
		pg->d[i][j]->Frame = pg->d[i][j]->Remainder;
		pg->d[i][j]->Remainder = f;
				
/* 		score=pg->d[i][j]->Donor->Score; */
/* 		class=pg->d[i][j]->Donor->class; */
/* 		pg->d[i][j]->Donor->Score = pg->d[i][j]->Acceptor->Score; */
/* 		pg->d[i][j]->Donor->class = pg->d[i][j]->Acceptor->class; */
/* 		pg->d[i][j]->Acceptor->Score = score; */
/* 		pg->d[i][j]->Acceptor->class = class; */

		/* Mark exon */
		pg->d[i][j]->selected = 0;
	      }   
	  }
    }
 
}

/* Restore frame/remainder in reverse exons read from gff file */
void UndoFrames(exonGFF* e, long n)
{
  long i;
  int f;
/*   float score; */
/*   short class; */
  /* Undo frame/rmd change*/
  for (i=0;i<n;i++)  
    if ((e+i)->Strand == '-')
      {
	f=(e+i)->Frame;
	(e+i)->Frame=(e+i)->Remainder;
	(e+i)->Remainder=f;
/* 	score=(e+i)->Donor->Score; */
/* 	class=(e+i)->Donor->class; */
/* 	(e+i)->Donor->Score = (e+i)->Acceptor->Score; */
/* 	(e+i)->Donor->class = (e+i)->Acceptor->class; */
/* 	(e+i)->Acceptor->Score = score; */
/* 	(e+i)->Acceptor->class = class; */

      }
 
}
