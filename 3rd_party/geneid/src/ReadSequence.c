/*************************************************************************
*                                                                        *
*   Module: ReadSequence                                                 *
*                                                                        *
*   Reading input sequences of DNA in fasta format                       *
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

/*  $Id: ReadSequence.c,v 1.8 2011/01/13 11:06:16 talioto Exp $  */

#include "geneid.h"

extern int VRB;

/* Predicting the length of a sequence before loading it */
long analizeFile(char* SequenceFile)
{
  struct stat *buf;
  long size;

  if ((buf = (struct stat*) malloc(sizeof(struct stat))) == NULL)
    printError("Not enough memory: analizing input sequence file"); 

  if (stat(SequenceFile,buf) != 0)
    printError("Impossible to access sequence file");

  /* Size (bytes) of this file */
  size = (long)buf->st_size;
  
  if (size == 0)
    printError("Empty input sequence file!");

  free(buf);

  return(size);
}

/* Get the header of the first DNA sequence (Name) */
/* The sequence file is allowed to contain more than one fasta sequence */
int IniReadSequence(FILE* seqfile, char* line)
{
  int res;
  char sAux[MAXSTRING];
  char cAux;

  /* Fasta format: >string1_string2...stringn \n ctacgatcgacg... */
  /* Searching ">" */
  res = fscanf(seqfile,"%c",&cAux);

  if ((res == -1) || cAux != '>')
	printError("Problems reading locusname (>)");

  /* Get locus name */
  res = fscanf(seqfile,"%s",sAux);
  if (res==-1)
	printError("Problems reading locusname");
  else
    strcpy(line,sAux);

  /* Jumping to the first fasta line (skipping the \n) */
  res = fscanf(seqfile,"%c",&cAux);

  while((cAux != '\n') && (res == 1))
    {
	  res = fscanf(seqfile,"%c",&cAux);
      if (res==-1)
		printError("Problems reading locusname");
    }
  
  if (res == EOF)
	printError("Problems reading locusname: unexpected end");
  
  return(res);
}

/* Reading content of current DNA sequence and the header of next one */
int ReadSequence (FILE* seqfile, char* Sequence, char* nextLocus)
{
  long pos;
  int res;
  char mess[MAXSTRING];
  char cAux;
  
  /* 1. Reading the current fasta sequence */
  /* fasta format = "atcgata...atta\n" */
  pos = 0;
  res = fscanf(seqfile,"%s\n",Sequence);
  while((res != EOF) && (Sequence[pos] != '>'))
    { 
      /* chars read = previous reading + current line */
      pos = pos + strlen(Sequence + pos);
      res = fscanf(seqfile,"%s",Sequence + pos);
      
      if (VRB && !(pos % MESSAGE_FREQ))
		{
          sprintf(mess, "...%ld",pos);
          printReadingInfo(mess);
		}
    }
  /* Repeat until the beginning of the next sequence is reached ">" */

  /* 2. Starting next Sequence */
  if (res != EOF)
    {
	  /* Forget every position after "pos" */
	  Sequence[pos]='\0';

      /* Delete '>' from the last read line */
	  nextLocus[0] = '\0';
      strcpy(nextLocus,Sequence+pos+1);

      /* Jumping until \n of the first fasta line */
      res = fscanf(seqfile,"%c",&cAux);
	  while((cAux != '\n') && (res == 1))
		{
		  res = fscanf(seqfile,"%c",&cAux);
		  if (res==-1)
			printError("Problems reading locusname");
		}
      
      if (res == EOF)
		printError("Problems reading locusname: unexpected end");
    }

  return(res);
}

