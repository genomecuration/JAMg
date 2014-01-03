/*  Copyright (c) 2003 by Mihaela Pertea.
 
* false.c = computes how many false positives/negatives we obtain
             by setting some threshold
*/


#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>


#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

int NO_ACC,NO_DON,NO_FA_ACC,NO_FA_DON;
float Scorep[100000],Scoren[500000];

void compute(int, int, int, float *, int *, float *, int *, float *);

int main ( int argc, char * argv [])
{ 
   FILE  * Infile, * Outfile;
   int i;
   float thr, pn, pp;
   int fn,fp;
   char s[1000];

   if  (argc < 9)
       {
        fprintf (stderr, "USAGE:  %s <InpScoreExFile> <InpScoreInFile> <OutFalseExFile> <OutFalseInFile> no_acc no_facc no_don no_fdon\n",
                    argv [0]);
        exit (EXIT_FAILURE);
       }   

   /* dealing with acceptors */ 

   NO_ACC=atoi(argv[5]);
   NO_FA_ACC=atoi(argv[6]);
   NO_DON=atoi(argv[7]);
   NO_FA_DON=atoi(argv[8]);

   Infile = fopen (argv [1], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
        exit (EXIT_FAILURE);
       }

   fgets(s,1000,Infile);

   for(i=0;i<NO_ACC;i++)
     fscanf(Infile,"%f",&Scorep[i]);

   fgets(s,1000,Infile);
   fgets(s,1000,Infile);

   for(i=0;i<NO_FA_ACC;i++)
     fscanf(Infile,"%f",&Scoren[i]);
       
   fclose(Infile);

   Outfile = fopen (argv [3], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [3]);
        exit (EXIT_FAILURE);
       }

   fprintf(Outfile," Threshold     False negatives     False positives\n");

   for(i=0;i<NO_ACC;i++)
     { 
       compute(i,NO_ACC,NO_FA_ACC,&thr,&fn,&pn,&fp,&pp);
       fprintf(Outfile,"%10.6f      %4d ( %5.2f%% )        %4d ( %5.2f%% )\n",thr,fn,pn,fp,pp);
     }

   fclose(Outfile);

   /* dealing with donors */ 

   Infile = fopen (argv [2], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
        exit (EXIT_FAILURE);
       }

   fgets(s,1000,Infile);

   for(i=0;i<NO_DON;i++){
     fscanf(Infile,"%f",&Scorep[i]);
   }

   fgets(s,1000,Infile);
   fgets(s,1000,Infile);

   for(i=0;i<NO_FA_DON;i++) {
     fscanf(Infile,"%f",&Scoren[i]);


   }
       
   fclose(Infile);

   Outfile = fopen (argv [4], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [4]);
        exit (EXIT_FAILURE);
       }

   fprintf(Outfile," Threshold     False negatives     False positives\n");

   for(i=0;i<NO_DON;i++)
     { 
       compute(i,NO_DON,NO_FA_DON,&thr,&fn,&pn,&fp,&pp);
       fprintf(Outfile,"%10.6f      %4d ( %5.2f%% )        %4d ( %5.2f%% )\n",thr,fn,pn,fp,pp);
     }

   fclose(Outfile);
       
 }

void compute(int where, int noa, int nofa, float *thr, int *fn, float *pn, int *fp, float *pp)
  {
    int i;

    *thr=Scorep[where];
    *fn=where;
    *pn=where*100;
    *pn=(*pn)/noa;

    i=0;
    while((Scorep[where]>Scoren[i])&&(i<nofa))
      i++;

    *fp=nofa-i;
    *pp=(*fp)*100;
    *pp=(*pp)/nofa;
  }
    

