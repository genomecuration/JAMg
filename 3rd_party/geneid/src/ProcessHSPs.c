/*************************************************************************
*                                                                        *
*   Module: ProcessHSPs                                                  *
*                                                                        *
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

extern float MRM;
extern int UTR;
extern int SRP;
extern float NO_SCORE;


/* Projection of HSPs: save the maximum for each nucleotide */
/* Requirement: HSPs must sorted by Position1 */
void HSPScan(packExternalInformation* external,
			 packHSP* hsp, 
			 int Strand, 
			 long l1, long l2)
{
  short x;
  short frameStart, frameEnd;
  long i,j;
  float scoreHSP;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = l1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] || 
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] || 
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			}
		}
	}
  else
	{
	  /* HSPs in REVERSE strand */
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  /* HSPs in reverse strand are reverse-sorted by Position2 */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j < l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = hsp->sPairs[x][i]->Score / 
					(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1);
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos2; 
					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1;
					  j--)
					{
					  if (scoreHSP > external->sr[x][j-l1] ||
						  external->sr[x][j-l1] == NO_SCORE)
						external->sr[x][j-l1] = scoreHSP;
					}
				}
			}
		}
	}
}
/* Projection of HSPs: save the maximum for each nucleotide */
/* Requirement: HSPs must sorted by Position1 */
void ReadScan(packExternalInformation* external,
			 packHSP* hsp, 
			 int Strand, 
			 long l1, long l2)
{
  short x;
  short frameStart, frameEnd;
  long i,j;
  float scoreHSP;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->readcount[x][i] = 0.0;
		  
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 < l1; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos2 >= l1) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 < l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
				  /* / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = l1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos1 >= l1, Pos2 <= l2-OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 <= l2-OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; 
/* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  /* Save the projection of HSPs into the array: including negative HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: right end is out (Pos2 > l2) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos1 <= l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score; 
/* 				    /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <= l2;
					  j++)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			}
		}
	}
  else
	{
	  /* HSPs in REVERSE strand */
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;

	  /* For each frame and strand, preprocess homology information */
	  /* HSPs in reverse strand are reverse-sorted by Position2 */
	  for(x=frameStart; x < frameEnd; x++)
		{
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->sr[x][i] = NO_SCORE;
		  /* Reset arrays: NO_SCORE values */
		  for(i=0; i < l2-l1+1; i++)
			external->readcount[x][i] = 0.0;
		  
		  if (hsp != NULL)
			{
			  /* A. Skip HSPs out of this range: [l1,l2] */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 > l2; 
				   i++)
				;
			  
			  /* B. Partial HSPs in this fragment: left end is out (Pos1 <= l2) */
			  for (i = external->iSegments[x]; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos2 > l2; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j < l2;
					  j++)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			  
			  /* C. Complete HSPs in this fragment (Pos2 <= l2,Pos1 >= l1+OVERLAP) */
			  for (; 
				   i < hsp->nSegments[x] && hsp->sPairs[x][i]->Pos1 >= l1+OVERLAP; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  /(hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* For each position in the HSP update the array sr */
				  for(j = hsp->sPairs[x][i]->Pos1; 
					  j <= hsp->sPairs[x][i]->Pos2 && j <l2;
					  j++)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			  
			  /* Update partial counter: previous HSPs are useless for next split */
			  external->iSegments[x] = i; 	  
			  
			  /* D. Partial HSPs in this fragment: left end is out (Pos1 < l1) */
			  for (; 
				   i < hsp->nSegments[x] &&  hsp->sPairs[x][i]->Pos2 >= l1; 
				   i++)
				{
				  /* Score value */
				  scoreHSP = (RREADS/MRM) * hsp->sPairs[x][i]->Score;
/* 				  / (hsp->sPairs[x][i]->Pos2 - hsp->sPairs[x][i]->Pos1 + 1); */
				  
				  /* Update the array sr with some positions of current HSPs */
				  for(j = hsp->sPairs[x][i]->Pos2; 
					  j >= hsp->sPairs[x][i]->Pos1 && j >= l1;
					  j--)
					{
					  if (external->sr[x][j-l1] == NO_SCORE){
					    external->sr[x][j-l1] = scoreHSP;
					  }else{
					    external->sr[x][j-l1] = MIN(external->sr[x][j-l1]+scoreHSP,COV);
					  }
					  external->readcount[x][j-l1] = external->readcount[x][j-l1] + hsp->sPairs[x][i]->Score;
					}
				}
			}
		}
	}
}

/* Preprocessing of HSPs projections */
void HSPScan2(packExternalInformation* external,
			  packHSP* hsp, 
			  int Strand, 
			  long l1, long l2)
{
  short x;
  long i;
  float previousScore;
  float previousReadCount;
  short frameStart, frameEnd;
  
  if (Strand == FORWARD)
    {
	  frameStart = 0; 
	  frameEnd = FRAMES;
	}
  else
	{
	  frameStart = FRAMES; 
	  frameEnd = 2*FRAMES;
	}
  
  for(x=frameStart; x < frameEnd; x++)
    {
      previousScore = 0.0;
      previousReadCount = 0.0;
      /* Screening the whole sequence to accumulate the sum in every base */
      for (i=l1; i<=l2; i++)
		{
		  /* Accumulating step */
		  if (UTR){
		    external->sr[x][i-l1] = previousScore + log(external->sr[x][i-l1] + 1);
		    previousScore = external->sr[x][i-l1];
		    external->readcount[x][i-l1] = previousReadCount + external->readcount[x][i-l1];
		    previousReadCount = external->readcount[x][i-l1];
		  }else{
		    external->sr[x][i-l1] = previousScore + external->sr[x][i-l1];
		    previousScore = external->sr[x][i-l1];
		  }
		}
    }  
}


/* Management function to score and filter exons */
void ProcessHSPs(long l1,
                long l2,
                int Strand,
		packExternalInformation* external,
                packHSP* hsp
                )
{

  /* Fill in the temporary HSP arrays (pre-processing) */
  /* GENIS hack */
  if (SRP)
	{
	  if (UTR){
	    printMess("Preprocessing read information: step 1");
	    ReadScan(external,hsp,Strand,l1,l2);
	  }else{
	    printMess("Preprocessing homology information: step 1");
	    HSPScan(external,hsp,Strand,l1,l2);
	  }

	  printMess("Preprocessing homology information: step 2");
	  HSPScan2(external,hsp,Strand,l1,l2);
	}
}







