/*************************************************************************
*                                                                        *
*   Module: CorrectExon                                                  *
*                                                                        *
*   Recompute exon limits according to its type                          *
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

/*  $Id: CorrectExon.c,v 1.8 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Fixing positions in the input sequence for the input exon */
void CorrectExon(exonGFF *e)
{
  /* Correct positions of begin / end of the exon */
  /* Stop codon included into the exon */
  e->offset1 =
    ((!strcmp(e->Type,"Terminal") || !strcmp(e->Type,"Single"))
	 && e->Strand =='-') ? 
    -LENGTHCODON + COFFSET : +COFFSET;
  
  e->offset2 = 
    ((!strcmp(e->Type,"Terminal") || !strcmp(e->Type,"Single"))
	 && e->Strand =='+')? 
    LENGTHCODON + COFFSET : +COFFSET;

}
/* Fixing positions in the input sequence for the input exon */
void CorrectUTR(exonGFF *e)
{
  int startOffset = 0;
  int donorOffset = 0;
  if (!strcmp(e->Type,sUTRFIRSTHALF)||!strcmp(e->Type,sUTR5INTERNALHALF)){
    donorOffset= -1;
    e->offset1 = (e->Strand =='-') ? -donorOffset + COFFSET : +COFFSET;
    e->offset2 = (e->Strand =='+')? donorOffset + COFFSET : +COFFSET;
  }
  if (!strcmp(e->Type,sUTRTERMINALHALF)||!strcmp(e->Type,sUTR3INTERNALHALF)){
      startOffset= 1 + LENGTHCODON;
      e->offset1 = (e->Strand =='+') ? startOffset + COFFSET : +COFFSET;
      e->offset2 = (e->Strand =='-')? -startOffset + COFFSET : +COFFSET;
    }
 
  
  /* Correct positions of begin / end of the exon */


 
}
/* Fixing positions in the input sequence for the ORF */
void CorrectORF(exonGFF* e) 
{
  /* Correct positions of begin / end of the ORF */
  /* Stop codon included into the ORF */
  e->offset1 = (e->Strand =='+') ? 
    COFFSET + COFFSET : -LENGTHCODON + COFFSET;
  
  e->offset2 = (e->Strand =='+')? 
    LENGTHCODON + COFFSET : COFFSET - COFFSET;
}

