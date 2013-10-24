/*************************************************************************
*                                                                        *
*   Module: PeakEdgeScore                                                *
*                                                                        *
*   Signal prediction by using a Position Weighted Array                 *
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

float PeakEdgeScore(long Position, 
				   int Strand, 
				   packExternalInformation* external, 
		    long l1, long l2, int win)
{
  int index;
  short trueFrame;
  long relPos;
  float Score = 0;
  float factor = 10;
/*   int win = 6; */
/*   char mess[MAXSTRING]; */

  relPos = Position - l1 + COFFSET;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame blast definition: according to the sequence start */
  trueFrame = 0;
  trueFrame += index;
  
  /* Access the sr array to obtain the homology score for current score */
  if (((relPos -win) >=0) && ((relPos +win) < LENGTHSi) && ((Position+win)<l2)){
    Score = (factor*(
		     ((external->sr[trueFrame][relPos + win] - external->sr[trueFrame][relPos])/win) 
		     - ((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -win])/win)
		     )
	     );
/*     sprintf(mess,"relPos: %ld\nslope 2: %f\nslope 1: %f",relPos,((external->sr[trueFrame][relPos +win] - external->sr[trueFrame][relPos])/win),((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -win])/win)); */
/*     printMess(mess); */
/*     sprintf(mess,"peakEdgeScore: %f",(factor*(((external->sr[trueFrame][relPos +win] - external->sr[trueFrame][relPos])/win)  */
/* 					    - ((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -win])/win)))); */
/*     printMess(mess); */
  
  }
  
  return(Score);
}

int ClusterEdge(long Position, 
				   int Strand, 
				   packExternalInformation* external, 
		    long l1, long l2)
{
  int index;
  short trueFrame;
  long relPos;
  int Score = 0;
  
/*   int win = 6; */
/*   char mess[MAXSTRING]; */

  relPos = Position - l1 + COFFSET;

  if (Strand == FORWARD)
    index = 0; 
  else
    index = FRAMES;
  
  /* Frame blast definition: according to the sequence start */
  trueFrame = 0;
  trueFrame += index;
  
  /* Access the sr array to obtain the homology score for current score */
  if (((relPos -2) >=0) && ((relPos +1) < LENGTHSi) && ((Position+1)<l2)){
    if (((external->sr[trueFrame][relPos -1] - external->sr[trueFrame][relPos -2]) == 0)&&((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -1]) > 0)){
      Score = 1;
    }
    if (((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -1]) > 0)&&((external->sr[trueFrame][relPos+1] - external->sr[trueFrame][relPos]) == 0)){
      Score = -1;
    }
    
/*     sprintf(mess,"relPos: %ld\nslope 2: %f\nslope 1: %f",relPos,((external->sr[trueFrame][relPos +win] - external->sr[trueFrame][relPos])/win),((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -win])/win)); */
/*     printMess(mess); */
/*     sprintf(mess,"peakEdgeScore: %f",(factor*(((external->sr[trueFrame][relPos +win] - external->sr[trueFrame][relPos])/win)  */
/* 					    - ((external->sr[trueFrame][relPos] - external->sr[trueFrame][relPos -win])/win)))); */
/*     printMess(mess); */
  
  }
  
  return(Score);
}
