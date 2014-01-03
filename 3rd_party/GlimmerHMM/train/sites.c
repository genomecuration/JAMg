//Copyright (c) 2003  by  Mihaela Pertea.

#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>

#define  TRUE  1
#define  FALSE  0
#define  ALPHABET_SIZE  4
#define  ACCEPTOR_LEN  29                    /* Positions +44,72 in a80 */
#define  ACCEPTOR_FILE_NAME "acc"
#define  ACCEPTOR_TREE_FILE "outex"
#define  ACCEPTOR_SIGNAL_OFFSET  24          /* Start of  AG  */
//#define  ACCEPTOR_THRESHOLD     -11.41341   //-6.84352

#define  ATG_LEN  19                    /* Positions +0,18 in a19 */
#define  ATG_SIGNAL_OFFSET  12          /* Start of  ATG  */
#define  ATG_FILE_NAME "atg.markov"
//#define  ATG_THRESHOLD -4.4 

#define  STOP_LEN  19                    /* Positions +0,18 in a19 */
#define  STOP_SIGNAL_OFFSET  4          /* Start of  ATG  */
#define  STOP_FILE_NAME "stop.markov"
//#define  ATG_THRESHOLD -4.4 
   
#define  DONOR_LEN  16                        /* Positions +5,20 in d80 */
#define  DONOR_FILE_NAME "don"
#define  DONOR_TREE_FILE "outin"
#define  DONOR_SIGNAL_OFFSET  5               /* Start of  GT  */
//#define  DONOR_THRESHOLD   -3.562148

#define  AD_LEN  40                        /* Positions +5,20 in d80 */
#define  AD_FILE_NAME  "ad.markov"
#define  AD_SIGNAL_OFFSET  18              /* Start of  GT  */
#define  AD_THRESHOLD  0.084694            /* For 4 false negatives */
   

#define  MARKOV_DEGREE  1
#define  MARKOV_LEN  4                     /* ALPHABET_SIZE ^ MARKOV_DEGREE */
#define  LOG_PROB_ACCEPTOR  log (1682.0 / 30191.0)
#define  LOG_PROB_NONACCEPTOR  log (28509.0 / 30191.0)
#define  LOG_PROB_DONOR  log (1682.0 / 30191.0)         /* Change this */
#define  LOG_PROB_NONDONOR  log (28509.0 / 30191.0)     /* Change this */
#define  LOW_SCORE  -99.0  /* Score if pattern does not have GT or AG signal */
#define  RETURN_TRUE_PROB  0

#define CODING_LEN 80

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

typedef unsigned int word;

extern char TRAIN_DIR[500];

int Is_AD(const int *, double *);
int  Is_Acceptor  (const int *, double *,double,int);
int  Is_Donor  (const int *, double *,double,int);
int  Is_Atg  (const int * , double *,double);
int  Is_Stop  (const int * , double *,double);
int  Acc  (const int *, double *,double, tree *t,int ind);
int  Don  (const int *, double *,double, tree *t,int ind);
int comp(const void *a, const void *b);
int findfile(const int * S, tree *t);
void readtree(char *line, tree *t, int start);
int find(char *line, int start);
//double func(double x);
//double prod(float *p1,double *p2,word dim);
int  Is_Cod_NonCod  (const int * , double *, int ind);

/*double prod(float *p1,double *p2,word dim)
{
  word i;
  double s=0;
  for(i=0;i<dim;i++) s+=((double)p1[i])*p2[i];
  return(s);
}

double func(double x)
{
 double f;
 f=1/(1+exp(-x));
 return(f);
}*/

//#define THR_ACC    -4.381725
//#define THR_ACC_EX -6.01535
//#define THR_ACC_IN -6.288323
//#define THR_DON    -3.084269
//#define THR_DON_EX -5.650019
//#define THR_DON_IN -6.85514

#define  Start_PosEx 56
#define  Stop_PosEx 84

#define  Start_PosIn 75
#define  Stop_PosIn 90

#define  Start_Cod 0
#define  Stop_Cod 79

#define Start_NoCod 82
#define Stop_NoCod 161


int Is_Acceptor(const int *B, double *Return_Score,double ACCEPTOR_THRESHOLD, int istacc)
{
  FILE  * Infile;
  double Score,S1,S2;
  static tree *tacc;
  static int readtacc=FALSE;
  char line[5000];
  int i,ind;
  int T[100];
  double score1,score2,score3;
  char accname[600];
  
  if(istacc && !readtacc) {
    
    
    // read the structure of the acceptor tree 

    sprintf(accname,"%s%s",TRAIN_DIR,ACCEPTOR_TREE_FILE);

    Infile = fopen (accname, "r");
    if  (Infile == NULL)
      {
	fprintf (stderr, "ERROR:  Unable to open file %s\n", ACCEPTOR_TREE_FILE);
	exit (EXIT_FAILURE);
      }
    
    tacc = (tree *) malloc(sizeof(tree));
    if (tacc == NULL) {fprintf(stderr," Memory allocation for tree failure.\n"); abort();}
    fgets(line, 5000, Infile);
    i=strlen(line);
    line[i-1]='\0';
    fclose(Infile);
    
    readtree(line, tacc, 0);
    readtacc=TRUE;
  }
   
  for(i=0;i<=Stop_PosEx-Start_PosEx;i++)
    T[i]=B[i+Start_PosEx];

  ind=Acc(T, &S1, ACCEPTOR_THRESHOLD, tacc,0);
  if(ind==0) return(0);
  if(istacc) Acc(T, &S2, ACCEPTOR_THRESHOLD, tacc,1);
  else S2=S1;
  score1=(S1+S2)/2;

  //  if(score1<=THR_ACC) score1=-99;
  
  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];

  Is_Cod_NonCod(T,&score2,0);

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];

  //  if(score2<=THR_ACC_EX) score2=-99;
  
  Is_Cod_NonCod(T,&score3,1);

  //  if(score3<=THR_ACC_IN) score3=-99;

  Score=score1+score2+score3;

  *Return_Score=Score;
  
  
  return Score >= ACCEPTOR_THRESHOLD;
	  
      
}  

int Is_Donor(const int *B, double *Return_Score, double DONOR_THRESHOLD, int istdon)
{
  FILE  * Infile;
  double Score,S1,S2;
  static tree *tdon;
  static int readtdon=FALSE;
  char line[5000];
  int ind,i;
  int T[100];
  double score1,score2,score3;
  char donname[600];

  if(istdon && !readtdon) {
    
   
    // read the structure of the donor tree 
    sprintf(donname,"%s%s",TRAIN_DIR,DONOR_TREE_FILE);
    Infile = fopen (donname, "r");
    if  (Infile == NULL)
      {
	fprintf (stderr, "ERROR:  Unable to open file %s\n", DONOR_TREE_FILE);
	exit (EXIT_FAILURE);
      }
    
    tdon = (tree *) malloc(sizeof(tree));
    if (tdon == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
    fgets(line, 5000, Infile);
    i=strlen(line);
    line[i-1]='\0';
    fclose(Infile);
    
    readtree(line, tdon, 0);
    readtdon=TRUE;
  }

  for(i=0;i<=Stop_PosIn-Start_PosIn;i++)
    T[i]=B[i+Start_PosIn];

  ind=Don(T, &S1, DONOR_THRESHOLD, tdon,0);
  if(ind==0) return(0);
  if(istdon) Don(T, &S2, DONOR_THRESHOLD, tdon,1);
  else S2=S1;
  score1=(S1+S2)/2;


  //  if(score1<=THR_DON) score1=-99;

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];

  Is_Cod_NonCod(T,&score2,2);

  //  if(score2<=THR_DON_EX) score2=-99;
  
  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];

  Is_Cod_NonCod(T,&score3,3);

  //  if(score3<=THR_DON_IN) score3=-99;

  Score=score1+score2+score3;

  *Return_Score=Score;
  
  return Score >= DONOR_THRESHOLD;
	  
      
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
   		if (t->left == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->left,n);
     }
	
 i+=2;n=i;
 if(line[i]=='(') 
 	{
 		i=find(line,i+1);
		t->right = (tree *) malloc(sizeof(tree));
   		if (t->right == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
        readtree(line,t->right,n);
     }
}

int find(char *line, int start)
{
 int stop,i;

 i=start;

 while(line[i]!=')')
 	if(line[i]=='(') i=find(line,i+1);
 	else i++;
 stop=i+1;
 return(stop);
}
 	

int comp(const void *a, const void *b)
{ 
  if(*(double *)a > *(double *)b) return(1);
  else if (*(double *)a==*(double *)b) return(0);
  else return(-1);

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


int  Acc  (const int * S, double * Return_Score, double ACCEPTOR_THRESHOLD, tree *t,int ind)

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
   char accname[600];
#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub, no;

/* see which acceptor you should use */

if(ind) 
  {
	no=findfile(S,t);
	sprintf(accname,"%s%s%d",TRAIN_DIR,ACCEPTOR_FILE_NAME,no);

  }
else 
  {
    strcpy(accname,TRAIN_DIR);
    strcat(accname,"acc1.mar");
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

   return(1);
  }



int  Don  (const int * S, double * Return_Score, double DONOR_THRESHOLD, tree *t,int ind)

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
   char donname[600];
   int no;

#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   /* see which donor file you should use */
if(ind)
   { no=findfile(S,t);
   sprintf(donname,"%s%s%d",TRAIN_DIR,DONOR_FILE_NAME,no);
   }
else 
{
  strcpy(donname,TRAIN_DIR);
  strcat(donname,"don1.mar");
  no=0;
}

   if  (! Tables_Loaded[no] )
       {

        Infile = fopen (donname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open donor file \"%s\"\n",
                        donname);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
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

        for  (i = MARKOV_DEGREE - 1;  i < DONOR_LEN;  i ++)
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

   return(1);
  }

  

int Is_AD(const int *S, double *Return_Score)
{
  FILE *Infile;
  static float Positive_Table[AD_LEN][ALPHABET_SIZE][MARKOV_LEN];
  static float Negative_Table[AD_LEN][ALPHABET_SIZE][MARKOV_LEN];
  static int Tables_Loaded=FALSE;
  double Positive_Sum, Negative_Sum, Score;
  char adname[600];

#if RETURN_TRUE_PROB
  double X,Y;
#endif
  int i,j,k,Ct,Sub;

  sprintf(adname,"%s%s",TRAIN_DIR,AD_FILE_NAME);

  if(!Tables_Loaded)
    {
      Infile=fopen(adname,"r");
      if(Infile == NULL)
	{
	  fprintf(stderr,"ERROR: Unable to open ads file \"%s\"\n",
		  adname);
	  exit(EXIT_FAILURE);
	}

      for(i=MARKOV_DEGREE-1;i<AD_LEN;i++)
	for(k=0;k<MARKOV_LEN;k++)
	  for(j=0;j<ALPHABET_SIZE;j++)
	    {
	      Ct=fscanf(Infile,"%f",&Positive_Table[i][j][k]);
	      if(Ct!=1)
		{
		  fprintf(stderr,"ERROR reading ads file \"%s\"\n",
                           adname);
                  exit(EXIT_FAILURE);
		}
	    }

      for(i=MARKOV_DEGREE-1;i<AD_LEN;i++)
	for(k=0;k<MARKOV_LEN;k++)
	  for(j=0;j<ALPHABET_SIZE;j++)
	    {
	      Ct=fscanf(Infile,"%f",&Negative_Table[i][j][k]);
	      if(Ct!=1)
		{
		  fprintf(stderr,"ERROR reading ad file \"%s\"\n",
			  adname);
		  exit(EXIT_FAILURE);
		}
	    }
      fclose(Infile);
      
      Tables_Loaded=TRUE;
    }
  
  if(S[AD_SIGNAL_OFFSET]!=0
     || S[AD_SIGNAL_OFFSET+1]!=2
     || S[AD_SIGNAL_OFFSET+2]!=2
     || S[AD_SIGNAL_OFFSET+3]!=3)
    {
      *Return_Score=LOW_SCORE;
      return FALSE;
    }

  Sub =0;
  for(i=0;i<MARKOV_DEGREE;i++)
    Sub=ALPHABET_SIZE*Sub+S[i];

  Positive_Sum=Positive_Table[MARKOV_DEGREE-1][0][Sub];
  Negative_Sum=Negative_Table[MARKOV_DEGREE-1][0][Sub];

  for(i=MARKOV_DEGREE;i<AD_LEN;i++)
    {
      j=S[i];
      Positive_Sum += Positive_Table[i][j][Sub];
      Negative_Sum += Negative_Table[i][j][Sub];
      Sub=ALPHABET_SIZE*(Sub%(MARKOV_LEN/ALPHABET_SIZE))+j;
    }

  Score=Positive_Sum-Negative_Sum;

#if RETURN_TRUE_PROB
  X=exp(Positive_Sum+LOG_PROB_ACCEPTOR);
  Y=exp(Negative_sum+LOG_PROB_NONACCEPTOR);
  *Return_Score=log(X/(X+Y));
#else
  *Return_Score=Score;
#endif
  
return Score >= AD_THRESHOLD;

}



int  Is_Atg  (const int * S, double * Return_Score, double ATG_THRESHOLD )

/* Evaluate string  S [0 .. (ATG_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely atg
*  site.  Also set  Return_Score  to the probability that it is an atg
*  site. */

{
   FILE  * Infile;
   static float  Positive_Table [ATG_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [ATG_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int  Tables_Loaded = FALSE;
   double  Positive_Sum, Negative_Sum, Score;
   int i, j, k, Ct, Sub;
   char atgname[600];

   if  (! Tables_Loaded)
       {
	sprintf(atgname,"%s%s",TRAIN_DIR,ATG_FILE_NAME); 
        Infile = fopen (atgname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open atg file \"%s\"\n",
                        ATG_FILE_NAME);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < ATG_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading atg file \"%s\"\n", 
                                atgname);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < ATG_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading atg file \"%s\"\n", 
                                ATG_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded = TRUE;
       }

   if  (S [ATG_SIGNAL_OFFSET] != 0
	   || S [ATG_SIGNAL_OFFSET + 1] != 3
           || S [ATG_SIGNAL_OFFSET + 2] != 2)    /* ATG */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < ATG_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [i] [j] [Sub];
      Negative_Sum += Negative_Table [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }

   Score = Positive_Sum - Negative_Sum;


   * Return_Score = Score;


   return  Score >= ATG_THRESHOLD;
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
   static  int  Tables_Loaded[4] = {FALSE,FALSE,FALSE,FALSE};
   double  Positive_Sum, Negative_Sum, Score, Threshold;
   char filename[600];
   int no;


#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Ct, Sub;

   no=ind;

   switch (no) {
   case 0: // case of exon in acceptor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"score_ex.acc");
     break;
   case 1: // case of intron in acceptor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"score_in.acc");
     break;
   case 2: // case of exon in donor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"score_ex.don");
     break;
   case 3: // case of intron in donor
     strcpy(filename,TRAIN_DIR);
     strcat(filename,"score_in.don");
     break;
   }

   if  (! Tables_Loaded[no] )
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

        Tables_Loaded [no] = TRUE;
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

   return (1);
  }




int  Is_Stop  (const int * S, double * Return_Score, double STOP_THRESHOLD )

/* Evaluate string  S [0 .. (STOP_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely stop
*  site.  Also set  Return_Score  to the probability that it is an stop
*  site. */

{
   FILE  * Infile;
   static float  Positive_Table [STOP_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static float  Negative_Table [STOP_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
   static int  Tables_Loaded = FALSE;
   double  Positive_Sum, Negative_Sum, Score;
   int i, j, k, Ct, Sub;
   char stopname[600];

   if  (! Tables_Loaded)
       {
	sprintf(stopname,"%s%s",TRAIN_DIR,STOP_FILE_NAME); 
        Infile = fopen (stopname, "r");
        if  (Infile == NULL)
            {
             fprintf (stderr, "ERROR:  Unable to open stop file \"%s\"\n",
                        STOP_FILE_NAME);
             exit (EXIT_FAILURE);
            }

        for  (i = MARKOV_DEGREE - 1;  i < STOP_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Positive_Table [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading stop file \"%s\"\n", 
                                stopname);
                    exit (EXIT_FAILURE);
                   }
              }

        for  (i = MARKOV_DEGREE - 1;  i < STOP_LEN;  i ++)
          for  (k = 0;  k < MARKOV_LEN;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Ct = fscanf (Infile, "%f", & Negative_Table [i] [j] [k]);
               if  (Ct != 1)
                   {
                    fprintf (stderr, "ERROR reading stop file \"%s\"\n", 
                                STOP_FILE_NAME);
                    exit (EXIT_FAILURE);
                   }
              }

        fclose (Infile);

        Tables_Loaded = TRUE;
       }

   if  ((S [STOP_SIGNAL_OFFSET] != 3 || S [STOP_SIGNAL_OFFSET + 1] != 0 || S [STOP_SIGNAL_OFFSET + 2] != 0) &&    /* TAA */
       (S [STOP_SIGNAL_OFFSET] != 3 || S [STOP_SIGNAL_OFFSET + 1] != 0 || S [STOP_SIGNAL_OFFSET + 2] != 2) &&    /* TAG */
       (S [STOP_SIGNAL_OFFSET] != 3 || S [STOP_SIGNAL_OFFSET + 1] != 2 || S [STOP_SIGNAL_OFFSET + 2] != 0))      /* TGA */
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }

   Sub = 0;
   for  (i = 0;  i < MARKOV_DEGREE;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Positive_Table [MARKOV_DEGREE - 1] [0] [Sub];
   Negative_Sum = Negative_Table [MARKOV_DEGREE - 1] [0] [Sub];

   for  (i = MARKOV_DEGREE;  i < STOP_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Positive_Table [i] [j] [Sub];
      Negative_Sum += Negative_Table [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
     }

   Score = Positive_Sum - Negative_Sum;


   * Return_Score = Score;


   return  Score >= STOP_THRESHOLD;
  }
