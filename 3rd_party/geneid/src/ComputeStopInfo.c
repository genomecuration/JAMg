/*************************************************************************
*                                                                        *
*   Module: ComputeStopInfo                                              *
*                                                                        *
*   Store information to avoid prediction of stops in coding frame       *
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

/* $Id: ComputeStopInfo.c,v 1.4 2011/01/13 11:06:16 talioto Exp $ */

#include "geneid.h"

/* Prevention of stop codons when assembling two consecutive exons */
void ComputeStopInfo(exonGFF* e, char* s)
{
  /* A. left Value: frame */
  if (e->Frame == 0)
	e->lValue = 0;
  else
	{
	  if (e->Frame == 1)
		{
		  if (s[e->Acceptor->Position] == 'A')
			e->lValue = 2;
		  else
			{
			  if (s[e->Acceptor->Position] == 'G')
				e->lValue = 3;
			  else
				e->lValue = 0;
			}
		}
	  else
		{
		  if (e->Frame == 2)
			{
			  if ((s[e->Acceptor->Position] == 'A' && s[e->Acceptor->Position+1] == 'G') ||
				  (s[e->Acceptor->Position] == 'A' && s[e->Acceptor->Position+1] == 'A') ||
				  (s[e->Acceptor->Position] == 'G' && s[e->Acceptor->Position+1] == 'A'))
				e->lValue = 1;
			  else
				e->lValue = 0;
			}
		}
	}

  /* B. right Value: remainder */  
  if (e->Remainder == 0)
	e->rValue = 0;
  else
	{
	  /* geneid remainder equals complement(true remainder) -> here it is rmd = 1 */
	  if (e->Remainder == 2)
		{
		  if (s[e->Donor->Position] == 'T')
			e->rValue = 1;
		  else
			e->rValue = 0;
		}
	  else
		{
		  /* geneid remainder equals complement(true remainder)-> here it is rmd = 2 */
		  if (e->Remainder == 1)
			{ 
			  if ((s[e->Donor->Position-1] == 'T') && (s[e->Donor->Position] == 'G'))
				e->rValue = 2;
			  else
				{
				  if ((s[e->Donor->Position-1] == 'T') && (s[e->Donor->Position] == 'A'))
					e->rValue = 3;
				  else
					e->rValue = 0;
				}
			}
		}
	}
}
