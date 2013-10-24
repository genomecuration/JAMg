/*************************************************************************
*                                                                        *
*   Module: SearchEvidenceExons                                          *
*                                                                        *
*   Extracting evidence exons between (l1, l2) for the current split     *
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

/*  $Id: SearchEvidenceExons.c,v 1.10 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Looking for annotations to include with current ab initio predictions */
void SearchEvidenceExons(packExternalInformation* p, 
						 packEvidence* evidence,
						 long l2)
{
  long i;
  
  /* Evidences with acceptor higher than l2 will be ignored */
  i = p->i1vExons;
  while ((i < evidence->nvExons) && ((evidence->vExons+i)->Acceptor->Position + (evidence->vExons+i)->offset1) <= l2)
    i++;
  
  p->i2vExons = i;
  
  p->ivExons = p->i2vExons - p->i1vExons;
}

/* Processing next block of annotations */
void SwitchCounters(packExternalInformation* p)
{
  p->i1vExons = p->i2vExons;
}

/* Reset counters in packEvidence for the next input sequence */
void resetEvidenceCounters(packExternalInformation* p)
{
  p->i1vExons = 0;
  p->i2vExons = 0;
}



