/*************************************************************************
*                                                                        *
*   Module: BuildU12Acceptors.c                                          *
*                                                                        *
*   Signal prediction by using a Position Weighted Array                 *
*   using U12 branch point                                               *
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

/* Function TRANS: char -> integer such that A=0, C=1, G=2 and T/U=2 */
extern int TRANS[];

/* Maximum allowed number of generic sites */
extern long NUMSITES;

/* Additional profiles */
extern int BP;
extern int PPT;
extern int UTR;

float ComputeU12BranchProfile(char* s,
			      long positionAcc,
			      long limitRight,
			      profile* p,
			      site* splicesite)
{
  float maxScore;
  float score;
  int index;
  int Opt;
  long end;
  long i,j;
  maxScore = -INF;
  /*      char mess[MAXSTRING];  */
  i = MAX(p->order,positionAcc - p->acc_context);
  end = MIN(positionAcc - p->dist + p->dimension - p->offset,limitRight);
  Opt = MAX(0,positionAcc - p->opt_dist);
  for (;
       i + p->dimension <= end;
       i++)
    {
      /* Applying the additional profile */
      score=0.0;
      for (j=0;j < p->dimension;j++)
	{
	  /* i is the position inside the region */
	  /* 5 is used because there are A,C,G,T and N */
	  index = OligoToInt(s + i + j - p->order, p->order+1,5);
		  
	  if (index >= p->dimensionTrans)
	    score = score + -INFI;
	  else
	    score = score + p->transitionValues[j][index];
	}
	   
      score = score - p->penalty_factor * (((float)(abs(i + p->offset -Opt))/((float)(p->acc_context - p->offset - p->opt_dist)))*((float)(abs(i + p->offset
																			       -Opt))/((float)(p->acc_context - p->offset - p->opt_dist)))); 
      if (score >= maxScore){
	maxScore = score;
	splicesite->PositionBP = i + p->offset - positionAcc;
      }
    }
  return maxScore;
}

float ComputePPTProfile(char* s,
			long positionAcc,
			long limitRight,
			profile* p,
			site* splicesite)
{
  float maxScore;
  float score;
  int index;
  long end;
  long i,j;

    
  maxScore = -INF;

  i = MAX(p->order,positionAcc + splicesite->PositionBP + 1);
  i = MAX(i,positionAcc - p->dist - p->dimension);
  end = MIN(positionAcc,limitRight);
  for (;
       i + p->dimension <= end;
       i++)
    {
      /* Applying the additional profile */
      score=0.0;
      for (j=0;j < p->dimension;j++)
	{
	  /* i is the position inside the region */
	  /* 5 is used because there are A,C,G,T and N */
	  index = OligoToInt(s + i + j - p->order, p->order+1,5);
		  
	  if (index >= p->dimensionTrans)
	    score = score + -INFI;
	  else
	    score = score + p->transitionValues[j][index];
	}
	  
      if (score >= maxScore){
	maxScore = score;
	splicesite->PositionPPT = i + p->offset - positionAcc;
      }
    }
  
  /* Cutoff for BranchPoint and PPtracts are useless */
  /* if (maxScore < p->cutoff) */
  /* 	maxScore = 0.0; */

  return maxScore;
}

/* Search for acceptor splice sites, using additional profiles */
long  BuildU12Acceptors(char* s,
			short class,
			char* type,
			char* subtype,
			profile* u12_p,
			profile* u12bp,
			profile* ppt,
			site* st, 
			long l1, 
			long l2,
			long ns,
			long nsites,
			int Strand,
			packExternalInformation* external) 
{ 
  int i,j;
  char* sOriginal;
  float score;
  float scoreBP;
  float scorePPT;
  float scoreAcc;
  /*   long ns,is; */
  long is;
  long left,right;
  int index;
  float cutoff;

  /* Back-up the origin of the sequence */
  sOriginal = s;
  
  /* Calculate a cutoff for the original acceptor signal prediction - results in fewer branch point calculations */
  cutoff = u12_p->cutoff + (U12ACC_CUTOFF_FACTOR * u12_p->bfactor);


  /* 1. Searching sites between beginning of the sequence and p->offset */
  if (!l1)
    {
      for (is = 0; is < u12_p->offset && (ns<NUMSITES); is++)
	{ 	  
	  if (ns<NUMSITES){
		  	  
	    scorePPT = 0.0;
	    scoreBP = 0.0;
	    scoreAcc = 0.0;
	    score=0.0;
	    /* Applying part of the profile */
	    for (i=u12_p->offset-is, j=0; i < u12_p->dimension; i++,j++) 
	      {
		/* i is the position inside the region */
		index = OligoToInt(s+j, u12_p->order+1,5);

		if (index >= u12_p->dimensionTrans)
		  score = score + -INFI;
		else
		  score = score + u12_p->transitionValues[i][index];
	      }
	    if (score >= cutoff) 
	      {
				  
		/* Using additional profiles */
		scoreBP = ComputeU12BranchProfile(sOriginal,u12_p->offset-is,l2,u12bp,&st[ns]);
		if (scoreBP >=u12bp->cutoff) {
		if (PPT)
		  scorePPT = ComputePPTProfile(sOriginal,u12_p->offset-is,l2,ppt,&st[ns]);
				
		scoreAcc = score;
		score = score + scoreBP;
		score = u12_p->afactor + (u12_p->bfactor * score); 
		 
		if(UTR){
		  score = score + PeakEdgeScore(is + u12_p->order,Strand,external,l1,l2,6);
		}
		if (score >=  u12_p->cutoff) 
		  {
		    st[ns].Position = is + u12_p->order;
		    st[ns].ScoreBP = scoreBP;
		    st[ns].ScorePPT = scorePPT;
		    st[ns].ScoreAccProfile = scoreAcc;
		    st[ns].Score=score;
		    st[ns].class= class;
		    /* st[ns].PositionBP= 5; */
		    strcpy(st[ns].type,type);
		    strcpy(st[ns].subtype,subtype);
		    ns++;
		  }
		}
	      }
	  }
	}
    }
  
  
  
  /* 2. Normal processing: predicting using the whole profile */
  /* left and right are the true boundaries of prediction */
  left  = MAX(0+u12_p->order, l1 - u12_p->offset);
  right = l2 - u12_p->offset;
  s += left;
  is = 0;     
  /* Case A: Using Markov chain with order 0: PWM */
  if (u12_p->order == 0)
    {
      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension-1) && (is < right- left + 1) && (ns<NUMSITES))
	{ 	
	  if (ns<NUMSITES){

	    scorePPT = 0.0;
	    scoreBP = 0.0;
	    scoreAcc = 0.0;
	    /* is = 0..right */
	    score=0.0;
	    for (i=0;i<u12_p->dimension;i++)
	      {
		/* i is the position inside the region */
		index = TRANS[(int)(*(s + i))];
		if (index >= u12_p->dimensionTrans)
		  score = score + -INFI;
		else
		  score = score + u12_p->transitionValues[i][index];
	      }
	    if (score >= cutoff) 
	      {
		/* Using additional profiles */
		scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);  
		if (scoreBP >=u12bp->cutoff) {	  
		if (PPT)
		  scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
		scoreAcc = score;
		score = score + scoreBP;
		score = u12_p->afactor + (u12_p->bfactor * score);
		if(UTR){
		  score = score + PeakEdgeScore(left + is + u12_p->offset,Strand,external,l1,l2,6);
		}
		if (score >= u12_p->cutoff) 
		  {
		    st[ns].Position = left + is + u12_p->offset;
		    st[ns].ScoreBP = scoreBP;
		    st[ns].ScorePPT = scorePPT;
		    st[ns].ScoreAccProfile = scoreAcc;
		    st[ns].Score=score;
		    st[ns].class= class;
		    /* st[ns].PositionBP= 5; */
		    strcpy(st[ns].type,type);
		    strcpy(st[ns].subtype,subtype);
		    ns++;
		  }
		}
	      }
	  }			
	  is++;
	  s++;
	}
    }
	
  /* case B: Using Markov chain with order 1: dinucleotides */
  else if (u12_p->order == 1)
    {

      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension-1) && (is < right- left + 1) && (ns<NUMSITES))
	{ 		
	  if (ns<NUMSITES){
	    /*Do for U12GTAG*/
			  
	    scorePPT = 0.0;
	    scoreBP = 0.0;
	    scoreAcc = 0.0;
	    /* is = 0..right */
	    score=0.0;
	    for (i=0;i<u12_p->dimension;i++)
	      {
		/* i is the position inside the region */
		index = 5*TRANS[(int)(*(s + i -1))] + TRANS[(int)(*(s + i))];
		if (index >= u12_p->dimensionTrans)
		  score = score + -INFI;
		else
		  score = score + u12_p->transitionValues[i][index];
	      }
			  
	    if (score >= cutoff) 
	      {
		/* Using additional profiles */
		scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);
		if (scoreBP >=u12bp->cutoff) {	  
		if (PPT)
		  scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
		scoreAcc = score;
		score = score + scoreBP;
		score = u12_p->afactor + (u12_p->bfactor * score); 
		if(UTR){
		  score = score + PeakEdgeScore(left + is + u12_p->offset,Strand,external,l1,l2,6);
		}

		if (score >= u12_p->cutoff) 
		  {
		    st[ns].Position = left + is + u12_p->offset;
		    st[ns].ScoreBP = scoreBP;
		    st[ns].ScorePPT = scorePPT;
		    st[ns].ScoreAccProfile = scoreAcc;
		    st[ns].Score=score;
		    st[ns].class= class;
		    /* st[ns].PositionBP= 5; */
		    strcpy(st[ns].type,type);
		    strcpy(st[ns].subtype,subtype);
		    ns++;
		  }
		}
	      }
	  }
	  is++;
	  s++;
	}
    }
  /* case C: Using Markov chain with order > 1 */
  else
    {
      /* discovering splice sites with current profile */
      while (*(s+u12_p->dimension-1) && (is < right- left + 1) && (ns<NUMSITES))
	{ 
	  if (ns<NUMSITES){
			  
	    scorePPT = 0.0;
	    scoreBP = 0.0;
	    scoreAcc = 0.0;
	    /* is = 0..right */
	    score=0.0;
	    for (i=0;i<u12_p->dimension;i++)
	      {
		/* i is the position inside the region */
		/* 5 is used because there are A,C,G,T and N */
		index = OligoToInt(s + i - u12_p->order, u12_p->order+1,5);
		if (index >= u12_p->dimensionTrans)
		  score = score + -INFI;
		else
		  score = score + u12_p->transitionValues[i][index];
	      }
	    if (score >= cutoff) 
	      {
		/* Using additional profiles */
		scoreBP = ComputeU12BranchProfile(sOriginal,left + is + u12_p->offset,l2,u12bp,&st[ns]);
		if (scoreBP >=u12bp->cutoff) {
		if (PPT)
		  scorePPT = ComputePPTProfile(sOriginal,left + is + u12_p->offset,l2,ppt,&st[ns]);
				
		scoreAcc = score;
		score = score + scoreBP;
		score = u12_p->afactor + (u12_p->bfactor * score); 
		if(UTR){
		  score = score + PeakEdgeScore(left + is + u12_p->offset,Strand,external,l1,l2,6);
		}

		if (score >= u12_p->cutoff) 
		  {
		    st[ns].Position = left + is + u12_p->offset;
		    st[ns].ScoreBP = scoreBP;
		    st[ns].ScorePPT = scorePPT;
		    st[ns].ScoreAccProfile = scoreAcc;
		    st[ns].Score=score;
		    st[ns].class= class;
		    strcpy(st[ns].type,type);
		    strcpy(st[ns].subtype,subtype);
		    ns++;
		  }
		}
	      }
	  }	
	  is++;
	  s++;
	}
    }
  if (ns >= nsites)
    printError("Too many predicted sites: decrease RSITES parameter");
  
  return(ns);  

}

 
  
