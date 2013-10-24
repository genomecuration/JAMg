/*************************************************************************
*                                                                        *
*   Module: DumpHash                                                     *
*                                                                        *
*   Hash table of exons used during backup of genes                      *
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

/*  $Id: DumpHash.c,v 1.9 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Maximum memory space to backup exons between 2 fragments of sequence */
extern long MAXBACKUPEXONS;

long MAXDUMPHASH;

/* Initialize the hash table */
void resetDumpHash(dumpHash* h)
{
  long i;
  
  MAXDUMPHASH = (long) (MAXBACKUPEXONS / HASHFACTOR);
  
  for(i=0; i<MAXDUMPHASH; i++)
    h->T[i] = NULL;
  h->total = 0;
}

/* Hash Function:: String -> Integer between 0..MAXDUMPHASH-1 */
long fDump(exonGFF* E)
{
  long sum;
  long sum2;
/*   long sum3; */
/*   long sum4; */
  long total; 
  long i;
    sum = (E->Frame * E->Acceptor->Position + E->Donor->Position); 
  
  for(i=0, sum2=0; i < strlen(E->Type); i++)
    sum2 = (i+1) * E->Type[i] + sum2;
  
  total = (sum + sum2 * E->Strand) % MAXDUMPHASH;
  
/*   sum = ((E->Frame + E->Acceptor->class) * E->Acceptor->Position + E->Donor->Position * (E->Donor->class + 1));  */
  
/*   for(i=0, sum2=0; i < strlen(E->Type); i++) */
/*     sum2 = (i+1) * E->Type[i] + sum2; */
/*   for(i=0, sum3=0; i < strlen(E->Donor->subtype); i++) */
/*     sum3 = (i+1) * E->Type[i] + sum3; */
/*   for(i=0, sum4=0; i < strlen(E->Acceptor->subtype); i++) */
/*     sum4 = (i+1) * E->Type[i] + sum4;   */
/*   total = (sum + sum2 + sum3 + sum2 * E->Strand) % MAXDUMPHASH; */
  
  return(total);
}

/* Save the new exon into the hash table */
void setExonDumpHash(exonGFF* E, dumpHash* h)
{
  dumpNode* p;
  dumpNode* n;
  long i;

  /* Computing hash value */
  i = fDump(E);
  
  /* Allocate the new node */
  if ((n = (dumpNode *)malloc(sizeof(dumpNode))) == NULL)
    printError("Not enough memory: dumpster node");
  
  /* Filling in the node with exon features */
  n->acceptor = E->Acceptor->Position;
  n->donor = E->Donor->Position;
  n->frame = E->Frame;
  n->strand = E->Strand;
  n->dclass = E->Donor->class;
  n->aclass = E->Acceptor->class;
/*   strcpy(n->dsub, E->Donor->subtype); */
/*   strcpy(n->asub, E->Acceptor->subtype); */
  strcpy(n->type,E->Type);

  n->exon = E;
  
  if(h->T[i] == NULL)
    {
      n->next = NULL;
      h->T[i] = n;
    }
  else
    {
      /* There are more nodes in this position: Colission */
      p = h->T[i];
      /* Insert at the begining of the list */
      h->T[i] = n;
      n->next = p;
    }
  
  /* How many exons have we got? */
  h->total++;
}

/* Finding an exon. Returns the address or NULL pointer */ 
exonGFF* getExonDumpHash(exonGFF* E, dumpHash* h)
{
  long i;
  int found=0;
  exonGFF* exon;
  dumpNode *p;

  exon = NULL;
  i = fDump(E);
  
  /* Nobody lives here */
  if(h->T[i]==NULL)
    exon = NULL;
  else
    {
      /* There are more nodes in this position */
      p = h->T[i];
      /* Searching until the first position free */
      while( p != NULL && !found )
		{
		  /* Searching on the synonymous list the wanted exon */
		  if ((p->acceptor == E->Acceptor->Position) &&
			  (p->donor == E->Donor->Position) &&
			  (p->frame == E->Frame) &&
			  (p->strand == E->Strand) &&
		      (p->dclass == E->Donor->class) &&
		      (p->aclass == E->Acceptor->class) &&
/* 		      (!strcmp(p->dsub,E->Acceptor->subtype))  && */
/* 		      (!strcmp(p->asub,E->Donor->subtype))  && */
		      (!strcmp(p->type,E->Type))  )
			{
			  found = 1;
			  exon = p->exon;
			}
		  p = p->next;
		}
      if (!found)
		exon = NULL;
    }
  return(exon);
}

/* Free memory of hash nodes (sinonimous) */
void freeDumpNodes(dumpNode* node)
{
  if (node == NULL)
    ;
  else
    {
      freeDumpNodes(node->next);
      free(node);
    }
}

/* Free memory of the whole hash table (dumpster is ready to work) */
void cleanDumpHash(dumpHash *h)
{
  long i;

  /* free all the exons of the hash table */
  for(i=0; i<MAXDUMPHASH; i++)
    {
      freeDumpNodes(h->T[i]);
      h->T[i] = NULL;
    }
  
  h->total = 0; 
}


