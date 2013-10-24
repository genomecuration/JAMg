/*************************************************************************
*                                                                        *
*   Module: Dictionary                                                   *
*                                                                        *
*   Implementation of a look-up table by using hash tables               *
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

/*  $Id: Dictionary.c,v 1.9 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

/* Initializing the dictionary: hash table and counter of keys */
void resetDict(dict* d)
{
  int i;
  
  for(i=0; i<MAXENTRY; i++)
    d->T[i] = NULL;
  d->nextFree = 0;
}

/* Hash Function:: String -> Integer between 0..MAXENTRY-1 */
int f(char s[])
{
  int i,total;
  
  for(i=0, total=0; i<strlen(s); i++)
    total = (i+1)*s[i] + total;
  total = total % MAXENTRY;
  
  return(total);
}

/* Assign a number-key to the new word and store it */
int setkeyDict(dict* d, char s[])
{
  int key;
  node *p;
  node *n;
  int i;
  
  /* If this word exists at the dictionary don't insert */
  key = getkeyDict(d,s);
  if (key == NOTFOUND)
    {
      i = f(s);
	  
      /* Alloc the new word */
      if ((n = (node *)malloc(sizeof(node))) == NULL)
		printError("Not enough memory: dictionary word");
	  
      /* Filling the node */
      strcpy(n->s,s);
      n->key = d->nextFree++;
      if(d->T[i] == NULL)
		{
		  n->next = NULL;
		  d->T[i] = n;
		}
      else
		{
		  /* There are more nodes in this position: Colission */
		  p = d->T[i];
		  /* Insert at the begining of the list */
		  d->T[i] = n;
		  n->next = p;
		}
      key = n->key;
    }
  return(key);
}


/* Shows the dictionary */
void showDict(dict* d)
{
  int i;
  node *p;
  
  printf("Dictionary: \n\n");
  for(i=0 ; i < MAXENTRY ; i++)
    {
      if(d->T[i]!=NULL)     
		{
		  /* There are more nodes in this position */
		  p = d->T[i];
		  /* Searching the first position free */
		  while( p!= NULL )
			{
			  printf("%-20s | \t\t\t %d\n",p->s,p->key);
			  p = p->next;
			  
			}
		}
    }    
}
/* Returns the key for the word request; NOTFOUND is Not found */ 
int getkeyDict(dict* d, char s[])
{
  int i;
  int found=0;
  int key;
  node *p;
  /* showDict(d); */
  key = NOTFOUND;

  /* Computing hash function */
  i = f(s);
  
  /* Empty list means not found */
  if(d->T[i]==NULL)
    key = NOTFOUND; 
  else
    {
      /* There are more nodes in this position: run the list */
      p = d->T[i];
      /* Searching until the first position not used */
      while( p != NULL && !found )
		{
		  /* Same hash value: compare to see if it is the same string */
		  if(!strcmp(s,p->s))
			{
			  found = 1;
			  key = p->key;
			}
		  p = p->next;
		}
      if(!found)
		key = NOTFOUND;
    }
  return(key);
}

/* Free memory of hash nodes (sinonimous) */
void freeNodes(pnode node)
{
  if (node == NULL)
    ;
  else
    {
      freeNodes(node->next);
      free(node);
    }
}

/* Free memory of the whole dictionary */
void freeDict(dict* d)
{
  int i;
  
  /* free all of the words in the dictionary */
  for(i=0; i<MAXENTRY; i++)
	{
      freeNodes(d->T[i]);
      d->T[i] = NULL;
	}
  
  /* free the dictionary */
  free(d); 
}

/* Binding the amino acid (key) to the new codon (word) */
void setAADict(dict* d, char s[], char aA)
{
  node* p;
  node* n;
  int i;
  
  i = f(s);

  /* Allocating the new word */
  if ((n = (node *)malloc(sizeof(node))) == NULL)
	printError("Not enough memory: AA dictionary node");
  
  /* Filling the node */
  strcpy(n->s,s);
  n->key = aA;
  if(d->T[i] == NULL)
	{
	  n->next = NULL;
	  d->T[i] = n;
	}
  else
	{
	  /* There are more nodes in this position: Colission */
	  p = d->T[i];
	  /* Insert at the begining of the list */
	  d->T[i] = n;
	  n->next = p;
	}
}

/* Returns the amino acid for the input codon; 'X' is Not found */ 
char getAADict(dict* d, char s[])
{
  int i;
  int found=0;
  int aa;
  node *p;

  aa = UNKNOWNAA;
  
  i = f(s);
  if(d->T[i]==NULL)
    aa = UNKNOWNAA;
  else
    {
      /* There are more nodes in this position */
      p = d->T[i];
      /* Searching until the first position free */
      while( p != NULL && !found )
		{
		  if(!strcmp(s,p->s))
			{
			  found = 1;
			  aa = p->key;
			}
		  p = p->next;
		}
      if(!found)
		aa = UNKNOWNAA;
    }
  return(aa);
}

