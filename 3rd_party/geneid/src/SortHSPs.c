/*************************************************************************
*                                                                        *
*   Module: SortHSP                                                      *
*                                                                        *
*   Sort by position Pos1 the input HSP list                             *
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

#include "geneid.h"

int Split(HSP** hsps,int i, int j)
{
  long k;
  long x,c1,c2,m;
  long pivot;
  long pivotPos1;
  HSP* tmp;
  
  /* Choosing the random pivot: always the first one */
/*   k = hsps[(i+j)/2]->Pos1; */
  k = hsps[i]->Pos1;
  /* How many elements have less value than pivot (value "k") */
  for(m = 0, x=i+1; x <=j; x++)
    if (hsps[x]->Pos1 <= k)
	m++;
  
  /* This will be the right place for the pivot */
  pivot = m+i;

  tmp = hsps[i];
  hsps[i] = hsps[pivot];
  hsps[pivot] = tmp;

  pivotPos1 = hsps[pivot]->Pos1; 
  c1 = i;
  c2 = pivot+1;
  while (c1 < m+i && c2 <= j) 
    if (hsps[c1]->Pos1 <= pivotPos1)
	  c1++;
    else
      {
		tmp = hsps[c1];
		hsps[c1] = hsps[c2];
		hsps[c2] = tmp;
		
		c2++;
      }

  return(pivot);
}

void quickSort(HSP** hsps, int i, int j)
{
  long pivot;

  if (i == j+1)
    /* Nothing */;
  else
    {
      pivot = Split(hsps,i,j);
      
      quickSort(hsps,i,pivot-1);
      quickSort(hsps,pivot+1,j);
    }
}

int RSplit(HSP** hsps,int i, int j)
{
  long k;
  long x,c1,c2,m;
  long pivot;
  long pivotPos1;
  HSP* tmp;
  
  /* Choosing the random pivot: always the first one */
/*   k = hsps[(i+j)/2]->Pos1; */
  k = hsps[i]->Pos1;
  /* How many elements have a higher value than pivot (value "k") */
  for(m = 0, x=i+1; x <=j; x++)
    if (hsps[x]->Pos1 >= k)
	m++;
  
  /* This will be the right place for the pivot */
  pivot = m+i;

  tmp = hsps[i];
  hsps[i] = hsps[pivot];
  hsps[pivot] = tmp;

  pivotPos1 = hsps[pivot]->Pos1; 
  c1 = i;
  c2 = pivot+1;
  while (c1 < m+i && c2 <= j) 
    if (hsps[c1]->Pos1 >= pivotPos1)
	  c1++;
    else
      {
		tmp = hsps[c1];
		hsps[c1] = hsps[c2];
		hsps[c2] = tmp;
		
		c2++;
      }

  return(pivot);
}

void RquickSort(HSP** hsps, int i, int j)
{
  long pivot;

  if (i == j+1)
    /* Nothing */;
  else
    {
      pivot = RSplit(hsps,i,j);
      
      RquickSort(hsps,i,pivot-1);
      RquickSort(hsps,pivot+1,j);
    }
}

void SortHSPs(packHSP* p)
{
  int frame;
  char mess[MAXSTRING];

  for (frame=0; frame < FRAMES; frame++)
	{
	  sprintf(mess,"Quicksorting FWD HSPs in frame %d",frame);
	  printMess(mess);
	  
	  quickSort(p->sPairs[frame],
				0,p->nSegments[frame]-1);

	  sprintf(mess,"\t%ld HSPs successfully quicksorted",p->nSegments[frame]);
	  printMess(mess);
	}

  for (frame=FRAMES; frame < 2*FRAMES; frame++)
	{
	  sprintf(mess,"R-Quicksorting RVS HSPs in frame %d",frame);
	  printMess(mess);
	  
	  RquickSort(p->sPairs[frame],
				 0,p->nSegments[frame]-1);

	  sprintf(mess,"\t%ld HSPs successfully r-quicksorted",p->nSegments[frame]);
	  printMess(mess);
	}
}

