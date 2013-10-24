/*************************************************************************
*                                                                        *
*   Module: input                                                        *
*                                                                        *
*   Formatted input of SSgff  (GFF and fasta)                            *
*                                                                        *
*   This file is part of the geneid Distribution                         *
*                                                                        *
*     Copyright (C) 2000 - Genis PARRA FARRE                             *
*                          Enrique BLANCO GARCIA                         *
*                          Roderic GUIGO SERRA                           * 
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

#include "SSgff.h"

extern int VRB;

FILE* OpenFile(char *FileName) 
{
  FILE *file;
  char mess[MAXSTRING];
  
  sprintf(mess,"File %s cannot be open for read",FileName);
  
  if ((file=fopen(FileName, "r"))==NULL)
    printError(mess);
  
  return(file);
}


long ReadSequence (FILE* seqfile, char* Locus, char* Sequence)
{
  long pos;
  int res;
  char cAux;

  /* line has the locus of the current sequence */
  res = fscanf(seqfile,">%s",Locus);  
  
    /* Jumping until \n of the first fasta line */
  res = fscanf(seqfile,"%c",&cAux); 
  while(cAux != '\n') 
    res = fscanf(seqfile,"%c",&cAux); 

    /* fasta format = "atcgata...atta\n" */
  pos = 0;
  res = fscanf(seqfile,"%s\n",Sequence);
  
    /* Only one sequence can be read */
     while(res != EOF)  
       {  
		 pos = pos + strlen(Sequence + pos);
		 res = fscanf(seqfile,"%s\n",Sequence + pos);
		 
		 if ( !(pos % 10000000) )
		   fprintf(stderr,"...%ld bp\n",pos);
		 if (pos >= SEQLENGTH)
		   printError("Not enough memory: change SEQLENGTH parameter ");

       }

     /* End of sequence */
    pos = pos + strlen(Sequence + pos); 
    return(pos); 
}



long ReadExonsGFF (FILE *file, exonGFF *exons, pack_gene *genes, dict* d, long *numgen, long Length)
{
  char line[MAXLINE];
  char c;
  long lastPos1;
  char mess[MAXSTRING];
 
  char *line1;
  char *line2;
  char *line3;
  char *line4;
  char *line5;
  char *line6;
  char *line7;
  char *line8;
  char *line9;

  long i,g;
  
  /* Coments: line begins with # */
  /* gff format = "Name  Source  Type  Begin  End  Score  Strand  Frame group */
  i = 0;
  g = 0;

  lastPos1 = -INFI;
  /* Skip comment lines */
  while(fgets(line,MAXLINE,file)!=NULL)
    {

      if(line[0] == '#' || line[0] == '\n')
        {
          /* Skip this line */
          printMess("Skipping comment line");
        }
      else
        {
          /* For each line extract the features (GFF format) */
          line1 = (char *) strtok(line,"\t");
          line2 = (char *) strtok(NULL,"\t");
          line3 = (char *) strtok(NULL,"\t"); 
          line4 = (char *) strtok(NULL,"\t");
          line5 = (char *) strtok(NULL,"\t");
          line6 = (char *) strtok(NULL,"\t");
          line7 = (char *) strtok(NULL,"\t");
          line8 = (char *) strtok(NULL,"\t");
          line9 = (char *) strtok(NULL,"\n");

	  
          if (line1 == NULL || line2 == NULL || line3 == NULL ||
              line4 == NULL || line5 == NULL || line6 == NULL ||
              line7 == NULL || line8 == NULL || line9 == NULL)
            {
              sprintf(mess, "Bad format: Exon GFF %ld\n",i);	      
              printError(mess);
            }
	  
          /* Exon Strand */
          if (sscanf(line7,"%c",&(exons[i].Strand))!= 1)
            {
              sprintf(mess, "Bad format Strand: Exon %ld\n",i);
              printError(mess);
            }
	  
          /* 1/2. Sequence and Source not used */
          /* 3. Exon Type */
          if (sscanf(line3,"%s",exons[i].Type) != 1)
            {
              sprintf(mess, "Bad format Type: Exon %ld\n",i);
              printError(mess);
            }
	      
          /* 4. Left position */
          if (sscanf(line4,"%ld",&(exons[i].Position1)) != 1)
            {
              sprintf(mess, "Bad format Position1: Exon %ld\n",i);
              printError(mess);
            }

          /* 4.b. Filename sort by position */
		  if (exons[i].Position1 < lastPos1) 
            {
              sprintf(mess, "Bad position(not sorted): Exon %ld\n",i);
              printError(mess); 
            }
          lastPos1 = exons[i].Position1;

          /* 5. Right position */
          if (sscanf(line5,"%ld",&(exons[i].Position2)) != 1)
            {
              sprintf(mess, "Bad format Position2: Exon %ld\n",i);
              printError(mess);
            }

          /* 5.b Boundaries position */
		  if ((exons[i].Position1) < 0 || (exons[i].Position1) > Length ||
			  (exons[i].Position2) < 0 || (exons[i].Position2) > Length)
			{
              sprintf(mess, "Coordinates out of range: Exon %ld\n",i);
              printError(mess);
			}
		  if ((exons[i].Position1) > (exons[i].Position2)) 
			{
              sprintf(mess, "Bad format positions : Exon %ld\n",i);
              printError(mess);
			}

          /* 6. Score = '.' or float */
          if (sscanf(line6,"%lf",&(exons[i].Score)) != 1)
            {
              if ((sscanf(line6,"%c",&c)!= 1) || (c!='.'))
                {
                  sprintf(mess, "Bad format Score: Exon %ld\n",i);
                  printError(mess);
                }
              exons[i].Score = NOSCORE;
            }
          
          /* 7. Strand was done before */
        
          /* 8. Frame = '.' or integer */
          if (sscanf(line8,"%hd",&(exons[i].Frame)) != 1)
            {
              if ((sscanf(line8,"%c",&c)!= 1) || (c!='.'))
                {
                  sprintf(mess, "Bad format Frame: Exon %ld\n",i);
                  printError(mess);
                }
              exons[i].Frame = NOFRAME;
            }     
          /* 9. Group */
          if (sscanf(line9,"%s",exons[i].Group) != 1)
            {
              sprintf(mess, "Bad format Group: Exon %ld\n",i);
              printError(mess);
            }
	      
		  /* Filling Gene structure */

          /* Assign a integer to the group name */
		  g = getkeyDict (d,(exons+i)->Group);

		  /* If it does not exists assign a new integer to the gene name */
		  if (g == NOTFOUND) 
			{
			  setkeyDict(d,(exons+i)->Group);
			  g = getkeyDict (d,(exons+i)->Group);
			  /* Initializing a new gene */
			  (genes+g)->First = (exons+i);
			  (genes+g)->numExons = 0;
			  (genes+g)->Last = (exons+i);
			  /* Increment the number of genes */
			  (*numgen)++;
			}
		  else 
			{
			  /* Assigning the to the previous exon of
				 this gene (genes+g) the NextExon variable
				 pointing to the current exon (exons+i)*/
			  (genes+g)->Last->NextExon = exons+i;
			  /* Increment numExons and assign the new Last exon*/
			  (genes+g)->numExons++;
			  (genes+g)->Last = (exons+i);
			}
		  
		  //printf ("%d %d %d %d \n",*numgen, (genes+g)->numExons,
		  //   (genes+g)->First->Position1,(genes+g)->Last->Position1);
			
          /* Process features from current exon */
          if (i > NUMEXONS)
            printError("Too many exons: Change NUMEXONS parameter");          
          i++;
	  
        }
    }
  return(i);
}


