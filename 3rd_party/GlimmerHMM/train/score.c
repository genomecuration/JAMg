/* Copyright (c) 2003  by  Mihaela Pertea. 
* sitesk.c compute a score for splice sites based on karlin's paper */

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
#define  ACCEPTOR_LEN  29                    
#define  ACCEPTOR_SIGNAL_OFFSET  24          /* Start of  AG  */
#define  ACCEPTOR_FILE_NAME "acc"
#define  ACCEPTOR_TREE_FILE "outex"
#define  ACCEPTOR_THRESHOLD 0

#define  DONOR_LEN  16                        
#define  DONOR_SIGNAL_OFFSET  5               /* Start of  GT  */
#define  DONOR_FILE_NAME "don"
#define  DONOR_TREE_FILE "outin"
#define  DONOR_THRESHOLD 0



#define  MARKOV_DEGREE  1
#define  MARKOV_LEN  4                     /* ALPHABET_SIZE ^ MARKOV_DEGREE */
#define  LOW_SCORE  -99.0  /* Score if pattern does not have GT or AG signal */
#define  RETURN_TRUE_PROB  0

#define CODING_LEN 80

#define MemCheck(X,Name) if (X == NULL) {fprintf(stderr,"%s: Memory allocation failure.\n",Name); abort();}

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

#define  Start_PosEx 56
#define  Stop_PosEx 84

#define  Start_PosIn 75
#define  Stop_PosIn 90

#define  Start_Cod 0
#define  Stop_Cod 79

#define Start_NoCod 82
#define Stop_NoCod 161


//#define ACC_SITES
//#define DON_SITES

typedef struct tree {
	int val;
   	int consens;
   	int poz;
   	int no;
   	struct tree *left;
   	struct tree *right;
   } tree;

const int  MAX_STRING_LEN = 400;
const int  NUM_POSITIONS = 162;



int  Is_Acceptor  (const int *, double *, tree *t,int ind);
int  Is_Donor  (const int *, double *, tree *t,int ind);
int  Is_Cod_NonCod  (const int * , double *, int ind);
int *basetoint(char sequence[], long length);
int comp(const void *a, const void *b);
int findfile(const int * S, tree *t);
void readtree(char *line, tree *t, int start);
int find(char *line, int start);

int main ( int argc, char * argv [])
{ 
   FILE  * Infile, *Outfile;
   char  S [MAX_STRING_LEN], T [MAX_STRING_LEN], Name[MAX_STRING_LEN];
   int *B;
   long i, N;
   long String_Len, Poz;
   double Score[500000],S1,S2;
   double score1,score2,score3;
   tree *tdon, *tacc;
   char line[5000];
   int istacc, istdon, message;

   if  (argc < 10)
       {
        fprintf (stderr, "USAGE:  %s <InpExFile> <InpExFaFile> <InpInFile> <InpInFaFile> <OutExFile> <OutInFile> <tacc> <tdon> <message>\n",
                    argv [0]);
        exit (EXIT_FAILURE);
       }   

   istacc=atoi(argv[7]);
   istdon=atoi(argv[8]);
   message=atoi(argv[9]);

   //#ifdef ACC_SITES

   /* read the structure of the acceptor tree */

   if(istacc) {
     Infile = fopen (ACCEPTOR_TREE_FILE, "r");
     if  (Infile == NULL)
       {
	 fprintf (stderr, "ERROR:  Unable to open file %s\n", ACCEPTOR_TREE_FILE);
	 exit (EXIT_FAILURE);
       }
     
     tacc = (tree *) malloc(sizeof(tree));
     if (tacc == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
     fgets(line, 5000, Infile);
     i=strlen(line);
     line[i-1]='\0';
     fclose(Infile);
     
     readtree(line, tacc, 0);
   }

   /* dealing with true acceptors */

   String_Len = 1 + Stop_PosEx - Start_PosEx;

   Infile = fopen (argv [1], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
        exit (EXIT_FAILURE);
       }

   if(message) printf("Scores for true acceptors\n");

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
       sscanf (T, "A%ld %s %s",&Poz, S,Name);
       //sscanf (T, "%*s %s %s", S,Name);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_PosEx, String_Len);
      T [String_Len] = '\0';
      B = basetoint(T,String_Len);
      
      Is_Acceptor(B, &S1, tacc,0);
      if(istacc) { 
	Is_Acceptor(B, &S2, tacc,1);
	score1=(S1+S2)/2;
      }
      else score1=S1;

      strncpy (T, S + Start_NoCod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);
      
      score2=0;
      Is_Cod_NonCod(B,&score2,0);
      
      strncpy (T, S + Start_Cod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score3=0;
      Is_Cod_NonCod(B,&score3,1);

      if(message) printf("%s %ld %.6f %.6f %.6f %.6f\n",Name,Poz,score1,score2,score3,score1+score2+score3);

      Score[N++]=score1+score2+score3;
      free(B);
    }  

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   Outfile = fopen (argv [5], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [5]);
        exit (EXIT_FAILURE);
       }

   fprintf(Outfile,"Scores for true acceptors\n");

   for(i=0;i<N;i++){
     fprintf(Outfile, " %.6f\n", Score[i]);
     }


   /* dealing with false acceptors */
       
   String_Len = 1 + Stop_PosEx - Start_PosEx;

   Infile = fopen (argv [2], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [2]);
        exit (EXIT_FAILURE);
       }

   if(message) printf("Scores for false acceptors\n");

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
      sscanf (T, "FA%ld %s %s",&Poz, S,Name);
      //sscanf (T, "%*s %s %s", S,Name);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_PosEx, String_Len);
      T [String_Len] = '\0';
      B = basetoint(T,String_Len);
      
      Is_Acceptor(B, &S1, tacc,0);
      if(istacc) {
	Is_Acceptor(B, &S2, tacc,1);
	score1=(S1+S2)/2;
      }
      else score1=S1;

      strncpy (T, S + Start_NoCod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score2=0;
      Is_Cod_NonCod(B,&score2,0);

      strncpy (T, S + Start_Cod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score3=0;
      Is_Cod_NonCod(B,&score3,1);


      if(message) printf("%s %ld %.6f %.6f %.6f %.6f\n",Name,Poz,score1,score2,score3,score1+score2+score3);

      Score[N++]=score1+score2+score3;

      free(B);
    }

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   fprintf(Outfile,"Scores for false acceptors\n");

   for(i=0;i<N;i++){
     fprintf(Outfile, " %.6f\n", Score[i]);
   }

   fclose(Outfile);

   //#endif ACC_SITES

   //#ifdef DON_SITES

   /*read the structure of the donor tree */

   if(istdon) {
     Infile = fopen (DONOR_TREE_FILE, "r");
     if  (Infile == NULL)
       {
	 fprintf (stderr, "ERROR:  Unable to open file %s\n", DONOR_TREE_FILE);
	 exit (EXIT_FAILURE);
       }
     
     tdon = (tree *) malloc(sizeof(tree));
     if (tdon == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
     fgets(line, 5000, Infile);
     i=strlen(line);
     line[i-1]='\0';
     fclose(Infile);
     
     readtree(line, tdon, 0);
   }


   /* dealing with true donors */

   String_Len = 1 + Stop_PosIn - Start_PosIn;

   Infile = fopen (argv [3], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [3]);
        exit (EXIT_FAILURE);
       }

   if(message) printf("Scores for true donors\n");

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
      sscanf (T, "D%ld %s %s", &Poz,S,Name);
       //sscanf (T, "%*s %s %s", S,Name);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_PosIn, String_Len);
      T [String_Len] = '\0';
      B = basetoint(T,String_Len);

      Is_Donor(B, &S1,tdon,0);
      if(istdon) {
	Is_Donor(B, &S2,tdon,1);
	score1=(S1+S2)/2;
      }
      else score1=S1;

      strncpy (T, S + Start_Cod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score2=0;
      Is_Cod_NonCod(B,&score2,2);

      strncpy (T, S + Start_NoCod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);
      
      score3=0;
      Is_Cod_NonCod(B,&score3,3);
	  
      if(message) printf("%s %ld %.6f %.6f %.6f %.6f\n",Name, Poz,score1,score2,score3,score1+score2+score3);      

      Score[N++]=score1+score2+score3;
      free(B);
    }

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   Outfile = fopen (argv [6], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [6]);
        exit (EXIT_FAILURE);
       }

   fprintf(Outfile,"Scores for true donors\n");

   for(i=0;i<N;i++)
   fprintf(Outfile, " %.6f\n", Score[i]);

   /* dealing with false introns */
       
  
   String_Len = 1 + Stop_PosIn - Start_PosIn;

   Infile = fopen (argv [4], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [4]);
        exit (EXIT_FAILURE);
       }

   if(message) printf("Scores for false donors\n");

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
      sscanf (T, "FD%ld %s %s",&Poz, S,Name);
      //sscanf (T, "%*s %s %s", S,Name);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_PosIn, String_Len);
      T [String_Len] = '\0';
      B = basetoint(T,String_Len);
      
      Is_Donor(B, &S1,tdon,0);
      if(istdon) {
	Is_Donor(B, &S2,tdon,1);
	score1=(S1+S2)/2;
      }
      else score1=S1;

      strncpy (T, S + Start_Cod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score2=0;
      Is_Cod_NonCod(B,&score2,2);

      strncpy (T, S + Start_NoCod, CODING_LEN);
      T [CODING_LEN] = '\0';
      B = basetoint(T,CODING_LEN);

      score3=0;
      Is_Cod_NonCod(B,&score3,3);


      if(message) printf("%s %ld %.6f %.6f %.6f %.6f\n",Name,Poz,score1,score2,score3,score1+score2+score3);      

      Score[N++]=score1+score2+score3;
      //Score[N++]=score1;

      free(B);
    }

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   fprintf(Outfile,"Scores for false donors\n");

   for(i=0;i<N;i++)
   fprintf(Outfile, " %.6f\n", Score[i]);

   fclose(Outfile);

   //#endif DON_SITES

 }
    
void readtree(char *line, tree *t, int start)
{
 int len;
 int i,n;
 char part[10];
 len=strlen(line);

 i=start;
 while((line[i]=='(')||(line[i]==' ')) i++;
 n=i;
 while(line[i]!=' ')
 {
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->val=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->consens=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->poz=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->no=atoi(part);

 t->left=NULL;
 t->right=NULL;

 i+=2;n=i;
 if(line[i]=='(') 
 	{
 		i=find(line,i+1);
		t->left = (tree *) malloc(sizeof(tree));
   		if (t->left == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->left,n);
     }
	
 i+=2;n=i;
 if(line[i]=='(') 
 	{
 		i=find(line,i+1);
		t->right = (tree *) malloc(sizeof(tree));
   		if (t->right == NULL) {fprintf(stderr,"%s: Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->right,n);
     }
}

int find(char *line, int start)
{
 int stop,i;



 i=start;

 while(line[i]!=')') {
 	if(line[i]=='(') i=find(line,i+1);
 	else i++;
 }
 stop=i+1;

 return(stop);
}
 	

int comp(const void *a, const void *b)
{ 
  if(*(double *)a > *(double *)b) return(1);
  else if (*(double *)a==*(double *)b) return(0);
  else return(-1);

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
      //fprintf(stderr,"non-acgt character in string (position %ld)\n",i);
      //intarray[i] = -1;
      intarray[i] = 1;
      break;
    }
  }

  return intarray;
}

int findfile(const int * S, tree *t)
{
	int val, cons, poz;


	val=t->val;
	cons=t->consens;

	if( cons !=-1)
	{ 
		poz=t->poz;
	    if(S[poz]==cons)
	    	val=findfile(S,t->left);
	    else val=findfile(S, t->right);
	}


	return(val);
}


int  Is_Acceptor  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (ACCEPTOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely acceptor
*  site.  Also set  Return_Score  to the probability that it is an acceptor
*  site. */

  {
   FILE  * Infile;
   static float  Positive_Table[300][ACCEPTOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table[300][ACCEPTOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int Tables_Loaded[300]={FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE}; 
   double  Positive_Sum, Negative_Sum, Score;
   char accname[20];
#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub, no;

/* see which acceptor you should use */

if(ind) 
  {
	no=findfile(S,t);
	sprintf(accname,"%s%d",ACCEPTOR_FILE_NAME,no);
  }
else 
  {
    strcpy(accname,"acc1.mar");
    no=0;
  }

   if  (! Tables_Loaded[no])
       {
        Infile = fopen (accname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open acceptor file \"%s\"\n",
                        accname);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < ACCEPTOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no][i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading acceptor file \"%s\"\n", 
                                ACCEPTOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < ACCEPTOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no][i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading acceptor file \"%s\"\n", 
                                ACCEPTOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded[no]  = TRUE;
       }

   if  (S [ACCEPTOR_SIGNAL_OFFSET] != 0
           || S [ACCEPTOR_SIGNAL_OFFSET + 1] != 2)    /* AG */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no][MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no][MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < ACCEPTOR_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
  


   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_ACCEPTOR);
   Y = exp (Negative_Sum + LOG_PROB_NONACCEPTOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

   return  Score >= ACCEPTOR_THRESHOLD;
  }



int  Is_Donor  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (DONOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */

  {
   FILE  * Infile;
   static float  Positive_Table [300][DONOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [300][DONOR_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int Tables_Loaded[300]={FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE}; 
   double  Positive_Sum, Negative_Sum, Score;
   char donname[20];
   int no;

#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   /* see which donor file you should use */
if(ind)
   { no=findfile(S,t);
    sprintf(donname,"%s%d",DONOR_FILE_NAME,no);
   }
else 
{
  strcpy(donname,"don1.mar");
  no=0;
}

   if  (! Tables_Loaded[no] )
       {

        Infile = fopen (donname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open donor file \"%s\" where no=%d\n",
                        donname,no);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\" where no=%d\n", 
                                DONOR_FILE_NAME,no);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\" where no=%d\n", 
                                DONOR_FILE_NAME,no);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded [no] = TRUE;
       }


   if  (S [DONOR_SIGNAL_OFFSET] != 2
           || S [DONOR_SIGNAL_OFFSET + 1] != 3)    /* GT */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < DONOR_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
 
   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_DONOR);
   Y = exp (Negative_Sum + LOG_PROB_NONDONOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

	if(Score==-99) printf("look one\n");

   return  Score >= DONOR_THRESHOLD;
  }




int  Is_Cod_NonCod  (const int * S, double * Return_Score, int ind)

/* Evaluate string  S [0 .. (CODING_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */

  {
   FILE  * Infile;
   static float  Positive_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static  int  Tables_Load[4] = {FALSE,FALSE,FALSE,FALSE};
   double  Positive_Sum, Negative_Sum, Score, Threshold;
   char filename[20];
   int no;


#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   no=ind;

   switch (no) {
   case 0: // case of exon in acceptor
     strcpy(filename,"score_ex.acc");
     Threshold = 0;
     break;
   case 1: // case of intron in acceptor
     strcpy(filename,"score_in.acc");
     Threshold = 0;
     break;
   case 2: // case of exon in donor
     strcpy(filename,"score_ex.don");
     Threshold = 0;
     break;
   case 3: // case of intron in donor
     strcpy(filename,"score_in.don");
     Threshold = 0;
     break;
   }

   if  (! Tables_Load[no] )
       {
        Infile = fopen (filename, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open donor file \"%s\"\n",
                        filename);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [no] [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading donor file \"%s\"\n", 
                                DONOR_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Load [no] = TRUE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [no] [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < CODING_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }
 


   Score = Positive_Sum - Negative_Sum;

#if  RETURN_TRUE_PROB
   X = exp (Positive_Sum + LOG_PROB_DONOR);
   Y = exp (Negative_Sum + LOG_PROB_NONDONOR);
   * Return_Score = log (X / (X + Y));
#else
   * Return_Score = Score;
#endif

	if(Score==-99) printf("look one\n");

   return  Score >= Threshold;
  }










