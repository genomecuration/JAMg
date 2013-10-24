/*************************************************************************
*                                                                        *
*   Module: account                                                      *
*                                                                        *
*   Accounting: results and run-time stats                               *
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

/*  $Id: account.c,v 1.8 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Init accounting (data structure) */
account* InitAcc()
{
  account* m; 

  /* Request memory to allocate acc. data */
  m = (account*) RequestMemoryAccounting();

  printMess("Reset accounting data");

  /* Reset counters */
  m->starts = 0;
  m->stops = 0;
  m->acc = 0;
  m->don = 0;
  m->tss = 0;
  m->tes = 0;
  m->starts_r = 0;
  m->stops_r = 0;
  m->acc_r = 0;
  m->don_r = 0;
  m->tss_r = 0;
  m->tes_r = 0;
  m->first = 0;
  m->internal = 0;
  m->terminal = 0;
  m->single = 0;
  m->orf = 0;
  m->zle = 0;
  m->utr = 0;
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;
  m->orf_r = 0;
  m->zle_r = 0;
  m->utr_r = 0;

  m->totalExons = 0;

  /* Computing the time used by every module (benchmarking): NOT USED */
  m->tSites = 0;
  m->tExons = 0;
  m->tGenes = 0;
  m->tSort = 0;
  m->tScore = 0;
  m->tBackup = 0;
  
  /* Starting the count (running time) */
  /* Real time */
  (void) time(&m->tStart);
  /* CPU time */
  clock();


  return(m);
}

/* Updating total number of sites and exons so far found in the sequence */
void updateTotals(account *m,
                  packSites* allSites,
                  packSites* allSites_r,
                  packExons* allExons,
                  packExons* allExons_r)
{
  m->starts += allSites->nStartCodons;
  m->stops += allSites->nStopCodons;
  m->acc += allSites->nAcceptorSites;
  m->don += allSites->nDonorSites;
  m->tss += allSites->nTS;
  m->tes += allSites->nTE;

  m->starts_r += allSites_r->nStartCodons;
  m->stops_r += allSites_r->nStopCodons;
  m->acc_r += allSites_r->nAcceptorSites;
  m->don_r += allSites_r->nDonorSites;
  m->tss_r += allSites_r->nTS;
  m->tes_r += allSites_r->nTE;
     
  m->first += allExons->nInitialExons;
  m->internal += allExons->nInternalExons;
  m->terminal += allExons->nTerminalExons;
  m->single += allExons->nSingles;
  m->orf += allExons->nORFs;
  m->zle += allExons->nZeroLengthExons;
  m->utr += allExons->nUtrInitialExons;
  m->utr += allExons->nUtrInitialHalfExons;
  m->utr += allExons->nUtrInternalExons;
  m->utr += allExons->nUtr5InternalHalfExons;
  m->utr += allExons->nUtr3InternalHalfExons;
  m->utr += allExons->nUtrTerminalHalfExons;
  m->utr += allExons->nUtrTerminalExons;

  m->first_r += allExons_r->nInitialExons;
  m->internal_r += allExons_r->nInternalExons;
  m->terminal_r += allExons_r->nTerminalExons;
  m->single_r += allExons_r->nSingles;  
  m->orf_r += allExons_r->nORFs;  
  m->zle_r += allExons->nZeroLengthExons;
  m->utr_r += allExons->nUtrInitialExons;
  m->utr_r += allExons->nUtrInitialHalfExons;
  m->utr_r += allExons->nUtrInternalExons;
  m->utr_r += allExons->nUtr5InternalHalfExons;
  m->utr_r += allExons->nUtr3InternalHalfExons;
  m->utr_r += allExons->nUtrTerminalHalfExons;
  m->utr_r += allExons->nUtrTerminalExons;

  m->totalExons = 
    m->first + 
    m->first_r +
    m->internal + 
    m->internal_r +
    m->terminal +
    m->terminal_r +
    m->single +
    m->single_r +
    m->orf +
    m->orf_r +
    m->zle +
    m->zle_r +
    m->utr +
    m->utr_r;
}

/* Reset acc. counters for the next input sequence */
void cleanAcc(account* m)
{
  
  /* Reset */
  m->starts = 0;
  m->stops = 0;
  m->acc = 0;
  m->don = 0;
  m->starts_r = 0;
  m->stops_r = 0;
  m->acc_r = 0;
  m->don_r = 0;

  m->first = 0;
  m->internal = 0;
  m->terminal = 0;
  m->single = 0;
  m->orf = 0;
  m->first_r = 0;
  m->internal_r = 0;
  m->terminal_r = 0;
  m->single_r = 0;
  m->orf_r = 0;

  m->totalExons = 0;
}

