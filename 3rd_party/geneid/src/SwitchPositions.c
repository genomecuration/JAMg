/*************************************************************************
*                                                                        *
*   Module: SwitchPositions                                              *
*                                                                        *
*   Exchanging left and right signals in reverse sense exons             *
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

/*  $Id: SwitchPositions.c,v 1.8 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Move left and right signals to preserve the ordering: left < right */
/* Exons are then pieces with special properties to be assembled */
void SwitchPositions(packExons* allExons)
{
  long i;
  site *c;
 
  for (i=0;i<allExons->nInitialExons;i++) 
    {
      c = (allExons->InitialExons+i)->Acceptor; 
      (allExons->InitialExons+i)->Acceptor = 
	(allExons->InitialExons+i)->Donor;
      (allExons->InitialExons+i)->Donor = c;
    }

  for (i=0;i<allExons->nInternalExons;i++) 
    {
      c = (allExons->InternalExons+i)->Acceptor; 
      (allExons->InternalExons+i)->Acceptor = 
	(allExons->InternalExons+i)->Donor;
      (allExons->InternalExons+i)->Donor = c;
    }

  for (i=0;i<allExons->nZeroLengthExons;i++) 
    {
      c = (allExons->ZeroLengthExons+i)->Acceptor; 
      (allExons->ZeroLengthExons+i)->Acceptor = 
	(allExons->ZeroLengthExons+i)->Donor;
      (allExons->ZeroLengthExons+i)->Donor = c;
    }

  for (i=0;i<allExons->nTerminalExons;i++) 
    {
      c = (allExons->TerminalExons+i)->Acceptor; 
      (allExons->TerminalExons+i)->Acceptor = 
	(allExons->TerminalExons+i)->Donor;
      (allExons->TerminalExons+i)->Donor = c;
    }

  for (i=0;i<allExons->nSingles;i++) 
    {
      c = (allExons->Singles+i)->Acceptor; 
      (allExons->Singles+i)->Acceptor = 
	(allExons->Singles+i)->Donor;
      (allExons->Singles+i)->Donor = c;
    }

  for (i=0;i<allExons->nORFs;i++) 
    {
      c = (allExons->ORFs+i)->Acceptor; 
      (allExons->ORFs+i)->Acceptor = 
	(allExons->ORFs+i)->Donor;
      (allExons->ORFs+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtrInitialExons;i++) 
    {
      c = (allExons->UtrInitialExons+i)->Acceptor; 
      (allExons->UtrInitialExons+i)->Acceptor = 
	(allExons->UtrInitialExons+i)->Donor;
      (allExons->UtrInitialExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtrInitialHalfExons;i++) 
    {
      c = (allExons->UtrInitialHalfExons+i)->Acceptor; 
      (allExons->UtrInitialHalfExons+i)->Acceptor = 
	(allExons->UtrInitialHalfExons+i)->Donor;
      (allExons->UtrInitialHalfExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtrInternalExons;i++) 
    {
      c = (allExons->UtrInternalExons+i)->Acceptor; 
      (allExons->UtrInternalExons+i)->Acceptor = 
	(allExons->UtrInternalExons+i)->Donor;
      (allExons->UtrInternalExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtr5InternalHalfExons;i++) 
    {
      c = (allExons->Utr5InternalHalfExons+i)->Acceptor; 
      (allExons->Utr5InternalHalfExons+i)->Acceptor = 
	(allExons->Utr5InternalHalfExons+i)->Donor;
      (allExons->Utr5InternalHalfExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtr3InternalHalfExons;i++) 
    {
      c = (allExons->Utr3InternalHalfExons+i)->Acceptor; 
      (allExons->Utr3InternalHalfExons+i)->Acceptor = 
	(allExons->Utr3InternalHalfExons+i)->Donor;
      (allExons->Utr3InternalHalfExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtrTerminalHalfExons;i++) 
    {
      c = (allExons->UtrTerminalHalfExons+i)->Acceptor; 
      (allExons->UtrTerminalHalfExons+i)->Acceptor = 
	(allExons->UtrTerminalHalfExons+i)->Donor;
      (allExons->UtrTerminalHalfExons+i)->Donor = c;
    }
  for (i=0;i<allExons->nUtrTerminalExons;i++) 
    {
      c = (allExons->UtrTerminalExons+i)->Acceptor; 
      (allExons->UtrTerminalExons+i)->Acceptor = 
	(allExons->UtrTerminalExons+i)->Donor;
      (allExons->UtrTerminalExons+i)->Donor = c;
    }

}

   
