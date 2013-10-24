/*************************************************************************
*                                                                        *
*   Module: RecomputePositions                                           *
*                                                                        *
*   Predicted signals in RVS: translation into FWD sense positions       *
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

/*  $Id: RecomputePositions.c,v 1.7 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Normalising positions on the reverse (and complemented) sequence */
void RecomputePositions(packSites* allSites, long l)
{
  long i;
  
  for (i=0; i < allSites->nStartCodons; i++)
    (allSites->StartCodons+i)->Position = 
      l-(allSites->StartCodons+i)->Position-1;

  for (i=0; i < allSites->nAcceptorSites; i++)
    (allSites->AcceptorSites+i)->Position = 
      l-(allSites->AcceptorSites+i)->Position-1;

  for (i=0; i < allSites->nDonorSites; i++)
    (allSites->DonorSites+i)->Position = 
      l-(allSites->DonorSites+i)->Position-1;

  for (i=0; i < allSites->nStopCodons; i++)
    (allSites->StopCodons+i)->Position = 
      l-(allSites->StopCodons+i)->Position-1;

  for (i=0; i < allSites->nTS; i++)
    (allSites->TS+i)->Position = 
      l-(allSites->TS+i)->Position-1;

  for (i=0; i < allSites->nTE; i++)
    (allSites->TE+i)->Position = 
      l-(allSites->TE+i)->Position-1;
}

   
