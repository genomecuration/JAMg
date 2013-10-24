/*************************************************************************
*                                                                        *
*   Module: SortSites                                                    *
*                                                                        *
*   Sort by position all of predicted sites                              *
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

/* Maximum allowed number of sites  */
extern long NUMSITES;


/* Struct for a node (list): pointer to site and to next node */
struct siteitem 
{
  site* Site;
  struct siteitem* nexitem;
} SiteItem;

/* Set free a list of nodes */
void FreeSiteItems(struct siteitem *q)
{
  if (q == NULL)
    ;
  else
    {
      FreeSiteItems(q->nexitem);
      free(q);
    }
}

/* Insert an site in the selected list according to its left position */
void UpdateSiteList(struct siteitem** p, site* InputSite)
{
  /* Insert the new node at the end of the list */
  while (*p!=NULL)
    p=&((*p)->nexitem); 
  
  /* Allocating new node for this site */
  if ((*p= (struct siteitem *) malloc(sizeof(struct siteitem))) == NULL)
    printError("Not enough memory: new siteitem (sorting)");
  
  /* Updating information for the node */
  (*p)->Site=InputSite;
  (*p)->nexitem=NULL;
}

/* Sort all of predicted sites by placing them in an array of lists */
/* corresponding every list to a beginning position for predicted sites */
void SortSites(site* Sites, long numSites,  site* sortedSites,     
               long l1, long l2)
{ 
  struct siteitem **SiteList, *q;
  long i;
  long pos;
  long n;
/*   int offset; */
  long l;
  long left;
  long right;
/*   long room; */
/*   char mess[MAXSTRING]; */
  long HowMany;

  /* 0. Creating the array for sorting: 1 - Length of fragment */
  left = l1;
  right = l2;
  l =  right - left + 1;

  /* Allocating memory for array SiteList */
  if ((SiteList = 
       (struct siteitem **) calloc(l + COFFSET + 10, sizeof(struct siteitem *))) == NULL)
    printError("Not enough memory: SiteList array (sorting)");

  /* Reset the positions, pointing to NULL */
  for (i=0;i<l + COFFSET;i++)
    SiteList[i]=NULL;

  /* 1. Insert sites in the proper list according to its left position */
  /* Adding predicted sites in the forward sense */
  for (i=0; i<numSites; i++) 
    {

      
      /* Correction between real and relative to fragment position */
      pos=(Sites+i)->Position - left + COFFSET;
      /* sprintf (mess,"site index: %ld abs pos: %ld rel pos, %ld", i, (Sites+i)->Position, pos); */
/*       printRes(mess); */
      /* Insert site in the proper list depending on the left position */
      UpdateSiteList(&(SiteList[pos]), Sites+i);
    }

  n = 0;

  /* 3. Traversing the table extracting the sites sorted by left position */
  HowMany = FSORT*NUMSITES;
  /* for (i=0, n=0; i<l; i++) */
  for (i=0; i<l+10; i++)
    {
      /* Extracting sites from list q */
      q=SiteList[i];
      while (q != NULL)
		{
		  /* Save the extracted site */

		  sortedSites[n].Position = q->Site->Position;
		  sortedSites[n].Score = q->Site->Score;
		  sortedSites[n].ScoreAccProfile = q->Site->ScoreAccProfile;
		  sortedSites[n].ScoreBP = q->Site->ScoreBP;
		  sortedSites[n].ScorePPT = q->Site->ScorePPT;
		  sortedSites[n].PositionBP = q->Site->PositionBP;
		  sortedSites[n].PositionPPT = q->Site->PositionPPT;
		  sortedSites[n].class = q->Site->class;
		  strcpy(sortedSites[n].subtype,q->Site->subtype);
		  strcpy(sortedSites[n].type,q->Site->type);
		  n++;
		  if (n >= HowMany)
			printError("Too many predicted sites: increase FSORT parameter");
	  
		  q=q->nexitem;
		}

	  /* Free chained items in the processed list q */
	  FreeSiteItems(SiteList[i]);       
	}

  /* Free empty array */
  free(SiteList);
  for (i=0; i<numSites; i++) 
    {  
      Sites[i].Position = sortedSites[i].Position;
      Sites[i].Score = sortedSites[i].Score ;
      Sites[i].ScoreAccProfile = sortedSites[i].ScoreAccProfile ;
      Sites[i].ScoreBP = sortedSites[i].ScoreBP ;
      Sites[i].ScorePPT = sortedSites[i].ScorePPT ;
      Sites[i].PositionBP = sortedSites[i].PositionBP ;
      Sites[i].PositionPPT = sortedSites[i].PositionPPT ;
      Sites[i].class = sortedSites[i].class ;
      strcpy(Sites[i].subtype,sortedSites[i].subtype );
      strcpy(Sites[i].type,sortedSites[i].type );		  
    }
}

   
