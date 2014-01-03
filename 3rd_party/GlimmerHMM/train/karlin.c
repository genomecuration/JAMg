// Copyright (c) 2003  by  Mihaela Pertea.

/* karlin.c train scores for splice sites based on karlin's paper */

#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>

#define  DEBUG  0
#define  TRUE  1
#define  FALSE  0
#define  ALPHABET_SIZE  4
#define  ACCEPTOR_LEN  23                    /* Positions +48,70 in a80 */
#define  ACCEPTOR_SIGNAL_OFFSET  20          /* Start of  AG  */
#define MemCheck(X,Name) if (X == NULL) {fprintf(stderr,"%s: Memory allocation failure.\n",Name); abort();}

#define  DONOR_LEN  9                        /* Positions +7,15 in d80 */
#define  DONOR_SIGNAL_OFFSET  3               /* Start of  GT  */

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

typedef struct tree {
	int val;
   	int consens;
   	int poz;
   	int no;
   	struct tree *left;
   	struct tree *right;
   } tree;

const int  MAX_STRING_LEN = 110;
const int  NUM_POSITIONS = 80;
int  Start_Position;
int  Stop_Position; 
int Nofile;
int off; /* the offset of gt or ag */

int C[100], V[100][4], M[100][100][4], W[110];

int *basetoint(char sequence[], long length);
void split(tree *t, char *Inname,char *Outname);
void save(FILE *Outfile, tree *t);
double chi_square(int cons, int poz);

int main ( int argc, char * argv [])
{ 
   tree *t;
   FILE  * Outfile;
   int i;

   if  (argc < 6)
       {
        fprintf (stderr, "USAGE:  %s <InpFile> <OutFileName> startpoz endpoz offset\n",
                    argv [0]);
        exit (EXIT_FAILURE);
       }

   /* off is the offset of the gt or ag in the DNA sequence */
   off=atoi(argv[5]);

   /* W  indicates for which position the consensus should be checked */
   for(i=0;i<MAX_STRING_LEN;i++)
   	W[i]=1;
   W[off]=0;
   W[off+1]=0;
      
   t = (tree *) malloc(sizeof(tree));
   if (t == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}

   Nofile=0;
   t->val=Nofile;
   t->consens=-1;
   t->poz=-1;
   t->no=-1;
   t->left=NULL;
   t->right=NULL;

   Start_Position = atoi(argv[3]);
   Stop_Position = atoi(argv[4]);

   split(t,argv[1],argv[2]);

   /* here I print the tree */
   Outfile = fopen (argv[2], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv[2]);
        exit (EXIT_FAILURE);
       }
   save(Outfile,t);
   fprintf(Outfile,"\n%d %d\n",Start_Position,Stop_Position);
   fclose(Outfile);

}

void save(FILE *Outfile, tree *t)
{
  fprintf(Outfile,"( %d %d %d %d l",t->val,t->consens, t->poz,t->no);
  if(t->left != NULL) save(Outfile, t->left);
  fprintf(Outfile," r");
  if(t->right != NULL) save(Outfile, t->right);
  fprintf(Outfile," )");  
}

void split(tree *t, char *Inname,char *Outname)
{
   FILE  *Infile, *NextInfile1, *NextInfile2, *Outfile;
   char  *S, *T ,*Copy;
   int *B;
   int i, j, k, N, N1, N2, maxi, tempnofile ;
   int String_Len;
   char Outname1[20], Outname2[20];
   double max, Sum;
   fpos_t pos;
   
   Infile = fopen (Inname, "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", Inname);
        exit (EXIT_FAILURE);
       }

   String_Len = 1 + Stop_Position - Start_Position;  

   /* initialize frequencies for consensus */
   for(i=0;i<String_Len;i++)
      	for(j=0;j<4;j++)
      		V[i][j]=0;

   T = (char *) malloc((MAX_STRING_LEN+1)*sizeof(char));
   MemCheck(T,"T");
   S = (char *) malloc((MAX_STRING_LEN+1)*sizeof(char));
   MemCheck(S,"S");

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
      sscanf (T, "%*s %s", S);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_Position, String_Len);
      T [String_Len] = '\0'; 
      N++;

      B = basetoint(T,String_Len); 

      /* compute frequencies for each position */
      
      for(i=0;i<String_Len;i++)
      	V[i][B[i]]++;
      	
      free(B); 
    }  

    t->no=N;

    /* compute consensus */
    for(i=0;i<String_Len;i++)
      {
      	max=0;
      	for(j=0;j<4;j++)
      		if(V[i][j]>max)
      		{
				C[i]=j;
				max=V[i][j];
			}
	  } 

	/*compute match between consensus and variable at position j */
	for(i=0;i<String_Len;i++)
		for(j=0;j<String_Len;j++)
			for(k=0;k<4;k++)
				M[i][j][k]=0;
	/*pos=0;
	  fsetpos(Infile,&pos);*/
	rewind(Infile);
	while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
      sscanf (T, "%*s %s", S);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_Position, String_Len);
      T [String_Len] = '\0'; 
      B = basetoint(T,String_Len);
      for(i=0;i<String_Len;i++)
      	if(C[i]==B[i])
      		for(j=0;j<String_Len;j++)
      				M[i][j][B[j]]++;
      free(B);
      }


    /* for each position compute the chi_square and where to do the partition */


   max=0;
   maxi=-1;
   for(i=0;i<String_Len;i++)
   if(W[i])
   {
	Sum=0;
	W[i]=0;
	for(j=0;j<String_Len;j++)
	{
		/*printf("i=%d,j=%d,Ci=%d,Vj0=%d,Vj1=%d,Vj2=%d,Vj3=%d,Mij0=%d,Mij1=%d,Mij2=%d,Mij3=%d\n",i,j,C[i],V[j][0],V[j][1],V[j][2],V[j][3],M[i][j][0],M[i][j][1],M[i][j][2],M[i][j][3]);
		fflush(stdout);*/

		if(W[j]) 
			Sum+=chi_square(i,j);
	}
	W[i]=1;

	if(Sum>max)
	{
		max=Sum;
		maxi=i;

	}
   } 

   printf("Consens %d pe poz %d\n",C[maxi],maxi);

   /* open output files */

   sprintf(Outname1,"%s%d",Outname,Nofile + 1);
   sprintf(Outname2,"%s%d",Outname,Nofile + 2);

   NextInfile1 = fopen (Outname1, "w");
   if  (NextInfile1 == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", Outname1);
        exit (EXIT_FAILURE);
       }

   NextInfile2 = fopen (Outname2, "w");
   if  (NextInfile2 == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", Outname2);
        exit (EXIT_FAILURE);
       }

   Copy = (char *) malloc((MAX_STRING_LEN+1)*sizeof(char));
   MemCheck(Copy,"Copy");
   N1=0;
   N2=0;
   /* split into 2 leaves if necessary */
   /*pos=0;
     fsetpos(Infile,&pos);*/
   rewind(Infile);
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {  
     	strcpy(Copy,T);
		sscanf (T, "%*s %s", S);
      	assert (strlen (S) == NUM_POSITIONS);
      	strncpy (T, S + Start_Position, String_Len);
      	T [String_Len] = '\0'; 
	    B = basetoint(T,String_Len); 
	    if(B[maxi] == C[maxi])
	    { 
	    	N1++;
	    	fprintf(NextInfile1,"%s",Copy);
	    }
	    else
	    {
			N2++;
	    	fprintf(NextInfile2,"%s",Copy);
	    }
	}

	free(S);
	free(T);
	free(Copy);

	fclose(NextInfile1);
	fclose(NextInfile2);

	printf("Frunze de %d si %d componente\n",N1,N2);

	tempnofile=Nofile;
	Nofile+=2;

	//if(N1+N2>224)
	//if(N1+N2>1000)
	if(N1>400 && N2>400)
	{   t->consens=C[maxi];
   		t->poz=maxi;   	
		t->left = (tree *) malloc(sizeof(tree));
   		if (t->left == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
		(t->left)->val=tempnofile+1;
   		(t->left)->consens=-1;
		(t->left)->poz=-1;
		(t->left)->no=N1;
   		(t->left)->left=NULL;
   		(t->left)->right=NULL;	
		if(N1>350) 
			{ 
				W[maxi]=0;
				split(t->left,Outname1,Outname);
			}
		W[maxi]=1;
		t->right = (tree *) malloc(sizeof(tree));
   		if (t->right == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
   		(t->right)->val=tempnofile+2;
   		(t->right)->consens=-1;
   		(t->right)->poz=-1;
   		(t->right)->no=N2;
   		(t->right)->left=NULL;
   		(t->right)->right=NULL;
		if(N2>350) split(t->right,Outname2,Outname);
	}
   		
   fclose(Infile);

}

/* compute X2(Ci,Xj) */
double chi_square(int c,int poz)
{
	double x[2],chi,cell1[2][4],cell2[2][4];
	int i,j,n;

	n=V[poz][0]+V[poz][1]+V[poz][2]+V[poz][3];
		
	for(j=0;j<4;j++)
	{
		cell1[0][j]=M[c][poz][j];
		cell1[1][j]=V[poz][j]-M[c][poz][j];
	}

	x[0]=M[c][poz][0]+M[c][poz][1]+M[c][poz][2]+M[c][poz][3];
	x[0]=x[0]/n;
	x[1]=(n-x[0]);
	x[1]=x[1]/n;

	for(i=0;i<2;i++)
		for(j=0;j<4;j++)
			cell2[i][j]=x[i]*V[poz][j];

	chi=0;

	for(i=0;i<2;i++)
		for(j=0;j<4;j++)
		{
			if(cell2[i][j]!=0)
			{ 
				chi+=(cell1[i][j]-cell2[i][j])*(cell1[i][j]-cell2[i][j])/cell2[i][j];
			}
		}

	return(chi);
}


/* convert the acgt sequence into a sequence of 0123 -- integers */
int *basetoint(char sequence[], long length)
{
  int *intarray;
  long i;
  
  intarray = (int *) malloc((length+1)*sizeof(int));
  MemCheck(intarray,"intarray");

  for(i = 0; i < length; i++) {
    switch(sequence[i]) {
    case 'a':
      intarray[i] = 0;
      break;
    case 'c':
      intarray[i] = 1;
      break;
    case 'g':
      intarray[i] = 2;
      break;
    case 't':
      intarray[i] = 3;
      break;
    default:
      fprintf(stderr,
	      "non-acgt character in string (position %ld)\n",
	      i);
      intarray[i] = -1;
      break;
    }
  }

  return intarray;
}
