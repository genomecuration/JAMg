/* Copyright (c) 2003  by  Mihaela Pertea. */

#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>

const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const int  MAX_LINE = 300;

#define  DEBUG  0
#define  TRUE  1
#define  FALSE  0
#define  ALPHABET_SIZE  4

#define  STOP_LEN  19                    /* Positions +0,18 in a19 */
#define  STOP_SIGNAL_OFFSET  4          /* Start of  STOP  */
#define  STOP_FILE_NAME "stop.markov"
#define  STOP_THRESHOLD -4.44 

#define  MARKOV_DEGREE  2
#define  MARKOV_LEN  16                     /* ALPHABET_SIZE ^ MARKOV_DEGREE */
#define  LOW_SCORE  -99.0  /* Score if pattern does not have GT or AG signal */
#define  RETURN_TRUE_PROB  0



#define MemCheck(X,Name) if (X == NULL) {fprintf(stderr,"%s: Memory allocation failure.\n",Name); abort();}

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

//#define THR
#define ACC_SITES
#define DON_SITES

typedef struct tree {
	int val;
   	int consens;
   	int poz;
   	int no;
   	struct tree *left;
   	struct tree *right;
   } tree;

const int  MAX_STRING_LEN = 400;
int  NUM_POSITIONS;
int MESSG;
int CODING_LEN;

int Start_Pos, Stop_Pos;

#define NO_OF_VAR 2 	/* number of variables for the decision trees */
#define NO_OF_TREES 10     /* number of decision trees  

/* this file contains freq's of all in-frame hexamers in all
   training sequences. */
#define IF_6MER_TRAIN        "exon.hexfreq"

/* this file contains the freq's of all hexamers in all training
   sequences, including exons, introns, and intergenic DNA */
#define TRAIN_6MERS          "all.hexfreq"

int no_of_trees = NO_OF_TREES;
int no_of_dimensions = NO_OF_VAR;
int no_of_categories = 2;
struct tree_node **cgtroots=NULL;
struct tree_node **gtnroots=NULL;
struct tree_node **nagroots=NULL;
struct tree_node **agcroots=NULL;
char TRAIN_DIR[500]="";

int  Is_Stop  (const int * , double *);

int *basetoint(char sequence[], long length);
int comp(const void *a, const void *b);

int main ( int argc, char * argv [])
{ 
   FILE  * Infile, *Outfile;
   char  S [MAX_STRING_LEN], T [MAX_STRING_LEN], Name[MAX_STRING_LEN];
   int *B;
   long i, N;
   long String_Len, Poz;
   double Score[500000];
   double score;
   char line[2000];
   double *rv;


   if  (argc < 5)
       {
        fprintf (stderr, "USAGE:  %s <InpFile> <InpFaFile> <OutFile> <NUM_POSITIONS> [<MESSG>]\n",
                    argv [0]);
        exit (EXIT_FAILURE);
       }   

   NUM_POSITIONS=atoi(argv[4]);
   if(argc==6) MESSG=atoi(argv[5]);
   else MESSG=0;

   Start_Pos = 0;
   Stop_Pos = NUM_POSITIONS-1;



   /* dealing with true stops */

   String_Len = 1 + Stop_Pos - Start_Pos;

   Infile = fopen (argv [1], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
        exit (EXIT_FAILURE);
       }

   if(MESSG) {
     printf("Scores for true stops\n");
   }


   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
       //      sscanf (T, "%*s %s %s", S,Name);
       sscanf (T, "A%ld %s %s",&Poz, S,Name);
      assert (strlen (S) == NUM_POSITIONS);
      strncpy (T, S + Start_Pos, String_Len);
      T [String_Len] = '\0';
      B = basetoint(T,String_Len);

      Is_Stop(B,&score);
      Score[N++]=score;

      if(MESSG) {
	printf("%s %ld %.6f\n",Name,Poz,score);
      }

      free(B);
    }  

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   Outfile = fopen (argv [3], "w");
   if  (Outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [3]);
        exit (EXIT_FAILURE);
       }

   fprintf(Outfile,"Scores for true stops\n");

   for(i=0;i<N;i++){
     fprintf(Outfile, " %.6f\n", Score[i]);
     }


   /* dealing with false acceptors */
       
   Infile = fopen (argv [2], "r");
   if  (Infile == NULL)
       {
        fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [2]);
        exit (EXIT_FAILURE);
       }

   if(MESSG) {
     printf("Scores for false stops\n");
   }

   N = 0;
   while  (fgets (T, MAX_STRING_LEN, Infile) != NULL)
     {
       //sscanf (T, "%*s %s %s", S,Name);
       sscanf (T, "FA%ld %s %s",&Poz, S,Name);
       assert (strlen (S) == NUM_POSITIONS);
       strncpy (T, S + Start_Pos, String_Len);
       T [String_Len] = '\0';
       B = basetoint(T,String_Len);
       
       Is_Stop(B,&score);
       Score[N++]=score;

       if(MESSG) {
	 printf("%s %ld %.6f\n",Name,Poz,score);
       }
       
       free(B);
     }

   fclose(Infile);

   qsort(Score, N, sizeof(double), comp);

   fprintf(Outfile,"Scores for false stops\n");

   for(i=0;i<N;i++){
     fprintf(Outfile, " %.6f\n", Score[i]);
   }

   fclose(Outfile);

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
      fprintf(stderr,
	      "non-acgt character in string (position %ld)\n",
	      i);
      intarray[i] = -1;
      break;
    }
  }

  return intarray;
}



void *  Safe_malloc  (size_t Len)

/* Allocate and return a pointer to  Len  bytes of memory.
*  Exit if fail. */

  {
   void  * P;

   P = malloc (Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  malloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_realloc  (void * Q, size_t Len)

/* Reallocate memory for  Q  to  Len  bytes and return a
*  pointer to the new memory.  Exit if fail. */

  {
   void  * P;

   P = realloc (Q, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  realloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }





char *  strdup  (char * & S1, const char * S2)

/* Allocate memory in  S1  for a copy of string  S2  and copy
*  it.  Return a pointer to  S1 . */

  {
   S1 = (char *) Safe_malloc (1 + strlen (S2));
   strcpy (S1, S2);

   return  S1;
  }



FILE *  File_Open  (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }



char  Complement  (char Ch)

/* Returns the DNA complement of  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      case  'r' :          // a or g
        return  'y';
      case  'y' :          // c or t
        return  'r';
      case  's' :          // c or g
        return  's';
      case  'w' :          // a or t
        return  'w';
      case  'm' :          // a or c
        return  'k';
      case  'k' :          // g or t
        return  'm';
      case  'b' :          // c, g or t
        return  'v';
      case  'd' :          // a, g or t
        return  'h';
      case  'h' :          // a, c or t
        return  'd';
      case  'v' :          // a, c or g
        return  'b';
      default :            // anything
        return  'n';
     }
  }



char  Filter  (char Ch)

//  Return a single  a, c, g or t  for  Ch .  Choice is to minimize likelihood
//  of a stop codon on the primary strand.

  {
   switch  (tolower (Ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  Ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        return  'c';
    }
  }



int  Read_String  (FILE * fp, char * & T, long int & Size, char Name [],
                   int Partial)

/* Read next string from  fp  (assuming FASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read.  If  Partial  is
*  true, then the first line must have 2 integers indicating positions
*  in the string and only those positions will be put into  T .  If
*  Partial  is false, the entire string is put into  T .  Sets  Name
*  to the first string after the starting '>' character. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct = FALSE;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch == EOF)
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }

int  Read_Multi_String  (FILE * fp, char * & T, long int & Size, char Name [],
                   int Partial, int no)

/* Read next string #no from  fp  (assuming MULTIFASTA format) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Partial indicates if first line has
*  numbers indicating a subrange of characters to read.  If  Partial  is
*  true, then the first line must have 2 integers indicating positions
*  in the string and only those positions will be put into  T .  If
*  Partial  is false, the entire string is put into  T .  Sets  Name
*  to the first string after the starting '>' character. */

  {
   char  * P, Line [MAX_LINE];
   long int  Len, Lo, Hi;
   int  Ch, Ct = FALSE;
   int i;

   i=0;

   while(i<no) {
     while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
       ;

     if  (Ch == EOF)
       return  FALSE;
     if (Ch == '>') i++;
   }

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   P = strtok (Line, " \t\n");
   if  (P != NULL)
       strcpy (Name, P);
     else
       Name [0] = '\0';
   Lo = 0;  Hi = LONG_MAX;
   if  (Partial)
       {
        P = strtok (NULL, " \t\n");
        if  (P != NULL)
            {
             Lo = strtol (P, NULL, 10);
             P = strtok (NULL, " \t\n");
             if  (P != NULL)
                 Hi = strtol (P, NULL, 10);
            }
        assert (Lo <= Hi);
       }

   Ct = 0;
   T [0] = '\0';
   Len = 1;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;
      if  (Ct < Lo || Ct > Hi)
          continue;

      if  (Len >= Size)
          {
           Size += INCR_SIZE;
           T = (char *) Safe_realloc (T, Size);
          }
      Ch = tolower (Ch);
      switch  (Ch)
        {
         case  'a' :
         case  'c' :
         case  'g' :
         case  't' :
         case  's' :
         case  'w' :
         case  'r' :
         case  'y' :
         case  'm' :
         case  'k' :
         case  'b' :
         case  'd' :
         case  'h' :
         case  'v' :
         case  'n' :
           break;
         default :
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, Name);
           Ch = 'n';
        }
      T [Len ++] = Ch;
     }

   T [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }


/************************************************************************/
/************************************************************************/




int  Is_Stop  (const int * S, double * Return_Score)

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


   if  (! Tables_Loaded)
       {
        Infile = fopen (STOP_FILE_NAME, "r");
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
                                STOP_FILE_NAME);
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

   if ((S [STOP_SIGNAL_OFFSET] != 3 || S [STOP_SIGNAL_OFFSET + 1] != 0 || S [STOP_SIGNAL_OFFSET + 2] != 0) &&    /* TAA */
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


