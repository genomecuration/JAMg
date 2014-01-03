/*  Copyright (c) 1999 by Mihaela Pertea, Arthur Delcher.


* Programmer:  A. Delcher
*     Written:  19 Jun 96
*    Modified:  16 Aug 96
*                2 Sep 96
*        File:  ~delcher/genes/build-markov-for.c
*
*  This program builds a markov model for the strings in the
*  specified input file moving in the forward (left-to-right)
*  direction.
*/


#include  <stdio.h>
#include  <math.h>
#include  <stdlib.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <time.h>
#include  <assert.h>


const int  DEBUG = 0;
const int  TRUE = 1;
const int  FALSE = 0;
const int  ALPHABET_SIZE = 4;
const int  MAX_STRING_LEN = 100000;
const int  NUM_POSITIONS = 80;
const int  MARKOV_DEGREE = 1;
const int  MARKOV_LEN = 4;             // ALPHABET_SIZE ^ MARKOV_DEGREE
const double  DEFAULT_COUNT = 0.0001;
const int  DEFAULT_START_POSITION = 0;
const int  DEFAULT_STOP_POSITION = 79;
int  Start_Position = 0;
int  Stop_Position = 79;
int  Append_Markov_Table = 0;
int num_pos;

#define  PURINE_EDIT_DIST  1

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif


int  Exact_Markov = 0;


void  Print_Consensus (float [NUM_POSITIONS] [ALPHABET_SIZE] [MARKOV_LEN], int);
char  Rev_Value  (int);
void *  Safe_malloc  (size_t);
void *  Safe_realloc  (void *, size_t);
void  Set_Counters (float [MARKOV_LEN], int);
int  Value  (char);


int main  (int argc, char * argv [])
{
  FILE  * Infile, * Outfile;
  char  S [MAX_STRING_LEN], T [MAX_STRING_LEN];
  char  * P;
  float  Table [NUM_POSITIONS] [ALPHABET_SIZE] [MARKOV_LEN], Sum;
  int  Ch_Ct [ALPHABET_SIZE] [NUM_POSITIONS];
  int  Sub;
  int  i, j, k, Ct, Max, N, String_Len;
  
  if  (argc < 3) {
    
    fprintf (stderr, "USAGE:  %s <input-file> <output-file>\n", argv [0]);
    exit (EXIT_FAILURE);
  }
  
  num_pos=NUM_POSITIONS;

  for  (i = 3;  i < argc;  i ++) {
   
    if (argv[i][0] == '_') {
      Ct = sscanf (argv [i] + 1, "%d", & num_pos);
      if  (Ct < 1) {
	
	fprintf (stderr,
		 "ERROR:  Illegal number of positions \"%s\"\n", argv [i]);
	exit (EXIT_FAILURE);
      }
      continue;
    }
    if  (argv [i] [0] == '+') {
      
      Ct = sscanf (argv [i] + 1, "%d,%d", & j, & k);
      if  (Ct < 1) {
	
	fprintf (stderr,
		 "ERROR:  Illegal range \"%s\"\n", argv [i]);
	exit (EXIT_FAILURE);
      }
      Start_Position = j;
      if  (Ct == 2)
	Stop_Position = k;
      if  (Start_Position < 0
	   || Stop_Position - Start_Position <= MARKOV_DEGREE
	   || Stop_Position >= num_pos) {
	
	fprintf (stderr,
		 "ERROR:  Illegal range \"%s\"\n", argv [i]);
	exit (EXIT_FAILURE);
      }
      continue;
    }
    if  (argv [i] [0] != '-') {
      
      fprintf (stderr, "ERROR:  Unrecognized parameter \"%s\"\n", argv [i]);
      exit (EXIT_FAILURE);
    }
    if  (strcmp (argv [i] + 1, "exact") == 0) {
    
      Exact_Markov = TRUE;
      printf ("Using exact model (no edit distance)\n");
    }
    else if  (strcmp (argv [i] + 1, "append") == 0) {
    
      Append_Markov_Table = TRUE;
    }
    else {
      
      fprintf (stderr, "ERROR:  Unrecognized option \"%s\"\n", argv [i]);
      exit (EXIT_FAILURE);
    }
  }

  String_Len = 1 + Stop_Position - Start_Position;

  for  (i = 0;  i < ALPHABET_SIZE;  i ++)
    for  (j = 0;  j < String_Len;  j ++)
      for  (k = 0;  k < MARKOV_LEN;  k ++)
	Table [j] [i] [k] = DEFAULT_COUNT;

  for  (i = 0;  i < ALPHABET_SIZE;  i ++)
    for  (j = 0;  j < String_Len;  j ++)
      Ch_Ct [i] [j] = 0;

  Infile = fopen (argv [1], "r");
  if  (Infile == NULL) {
    
    fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [1]);
    exit (EXIT_FAILURE);
  }
  if  (Append_Markov_Table)
    Outfile = fopen (argv [2], "a");
  else
    Outfile = fopen (argv [2], "w");
  if  (Outfile == NULL) {
    
    fprintf (stderr, "ERROR:  Unable to open file %s\n", argv [2]);
    exit (EXIT_FAILURE);
  }

  N = 0;
  while  (fgets (T, MAX_STRING_LEN, Infile) != NULL) {
  
    N ++;
    sscanf (T, "%*s %s", S);
    assert (strlen (S) == num_pos);
    strncpy (T, S + Start_Position, String_Len);
    T [String_Len] = '\0';

    Sub = 0;
    for  (i = 0;  i < MARKOV_DEGREE;  i ++)
      Sub = ALPHABET_SIZE * Sub + Value (T [i]);
    Set_Counters (Table [MARKOV_DEGREE - 1] [0], Sub);
    
    for  (i = MARKOV_DEGREE;  i < String_Len;  i ++) {
      
      j = Value (T [i]);
      Set_Counters (Table [i] [j], Sub);
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + j;
    }

    for  (i = 0;  i < String_Len;  i ++)
      Ch_Ct [Value (T [i])] [i] ++;
  }
   
  fclose (Infile);

  Sum = 0.0;
  for  (i = 0;  i < MARKOV_LEN;  i ++)
    Sum += Table [MARKOV_DEGREE - 1] [0] [i];
  for  (i = 0;  i < MARKOV_LEN;  i ++)
    Table [MARKOV_DEGREE - 1] [0] [i]
      = log (Table [MARKOV_DEGREE - 1] [0] [i] / Sum);

  for  (j = MARKOV_DEGREE;  j < String_Len;  j ++)
    for  (k = 0;  k < MARKOV_LEN;  k ++) {
       
      Sum = 0.0;
      for  (i = 0;  i < ALPHABET_SIZE;  i ++)
	Sum += Table [j] [i] [k];
      for  (i = 0;  i < ALPHABET_SIZE;  i ++)
	Table [j] [i] [k] = log (Table [j] [i] [k] / Sum);
    }

  for  (i = MARKOV_DEGREE - 1;  i < String_Len;  i ++)
    for  (k = 0;  k < MARKOV_LEN;  k ++) {
      for  (j = 0;  j < ALPHABET_SIZE;  j ++)
	fprintf (Outfile, " %.6f", Table [i] [j] [k]);
      fprintf (Outfile, "\n");
    }

  fclose (Outfile);

  printf ("Simple Consensus = ");
  for  (i = 0;  i < String_Len;  i ++) {
    
    k = 0;
    Max = Ch_Ct [0] [i];
    for  (j = 1;  j < ALPHABET_SIZE;  j ++)
      if  (Ch_Ct [j] [i] > Max) {
	
	k = j;
	Max = Ch_Ct [j] [i];
      }
    putchar (Rev_Value (k));
  }
  putchar ('\n');

  printf ("Markov Consensus = ");
  Print_Consensus (Table, String_Len);
  putchar ('\n');
  
  return  0;
}



void  Print_Consensus (float T [NUM_POSITIONS] [ALPHABET_SIZE] [MARKOV_LEN],
                       int String_Len)

/* Print the consensus sequence determined by  T  representing stings
*  of length  String_Len . */

  {
   int  i, j, k, q, Sub, End_Position;
   float  Max, X, M [MARKOV_LEN] [NUM_POSITIONS];
   int  From [MARKOV_LEN] [NUM_POSITIONS];
   char  String [81];
   
   String [String_Len] = '\0';

   for  (i = 0;  i < MARKOV_LEN;  i ++)
     M [i] [MARKOV_DEGREE - 1]
                   = T [MARKOV_DEGREE - 1] [0] [i];

   for  (i = MARKOV_DEGREE;  i < String_Len;  i ++)
     {
      for  (Sub = 0;  Sub < MARKOV_LEN;  Sub ++)
        {
         Max = M [Sub / ALPHABET_SIZE] [i - 1]
                    + T [i] [Sub % ALPHABET_SIZE] [Sub / ALPHABET_SIZE];
         k = Sub / ALPHABET_SIZE;
         for  (j = 1;  j < ALPHABET_SIZE;  j ++)
           {
            q = Sub / ALPHABET_SIZE + j * (MARKOV_LEN / ALPHABET_SIZE);
            X = M [q] [i - 1] + T [i] [Sub % ALPHABET_SIZE] [q];
            if  (X > Max)
                {
                 Max = X;
                 k = q;
                }
           }
         M [Sub] [i] = Max;
         From [Sub] [i] = k;
        }
     }

   End_Position = String_Len - 1;
   k = 0;
   Max = M [k] [End_Position];
   for  (j = 1;  j < MARKOV_LEN;  j ++)
     if  (Max < M [j] [End_Position])
         {
          Max = M [j] [End_Position];
          k = j;
         }
   for  (i = End_Position;  i >= MARKOV_DEGREE;  i --)
     {
      String [i] = Rev_Value (k % ALPHABET_SIZE);
      k = From [k] [i];
     }
   for  ( ;  i >= 0;  i --)
     {
      String [i] = Rev_Value (k % ALPHABET_SIZE);
      k /= ALPHABET_SIZE;
     }

   printf ("%s\n", String);


   Sub = 0;
   Max = T [MARKOV_DEGREE - 1] [0] [0];
   for  (i = 1;  i < MARKOV_LEN;  i ++)
     if  (T [MARKOV_DEGREE - 1] [0] [i] > Max)
         {
          Sub = i;
          Max = T [MARKOV_DEGREE - 1] [0] [i];
         }

   j = Sub;
   for  (i = MARKOV_DEGREE - 1;  i >= 0;  i --)
     {
      String [i] = Rev_Value (j % ALPHABET_SIZE);
      j /= ALPHABET_SIZE;
     }

   for  (i = MARKOV_DEGREE;  i < String_Len;  i ++)
     {
      k = 0;
      Max = T [i] [0] [Sub];
      for  (j = 1;  j < ALPHABET_SIZE;  j ++)
        if  (T [i] [j] [Sub] > Max)
            {
             k = j;
             Max = T [i] [j] [Sub];
            }

      String [i] = Rev_Value (k);
      Sub = ALPHABET_SIZE * (Sub % (MARKOV_LEN / ALPHABET_SIZE)) + k;
     }

   printf ("******** Old Way = %s", String);

   return;
  }



char  Rev_Value  (int n)

/* Returns the character represented by the value  n . */

  {
   switch  (n)
     {
      case  0 :
        return  'a';
      case  1 :
        return  'c';
      case  2 :
        return  'g';
      case  3 :
        return  't';
      default :
        return  '?';
     }
  }



void *  Safe_malloc  (size_t Len)

/* Same as  malloc  except prints error and halts if
*  unsuccessful. */

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

/* Same as  realloc  except prints error and halts if
*  unsuccessful. */

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



void  Set_Counters (float Ct [MARKOV_LEN], int Sub)

/* Add 1 to all counters in  Ct at subscripts equal to  Sub
*  or 1 base-ALPHABET_SIZE digit off from  Sub . */

  {
   int  i, j, A, Div_By, Mod_By;

   Ct [Sub] += 1.0;

   if  (! Exact_Markov)
       {
        j = Sub % ALPHABET_SIZE;
        A = Sub - j;
#if  PURINE_EDIT_DIST
//        Ct [A + (j + 2) % ALPHABET_SIZE] += 1.0;
#else
        for  (i = 0;  i < ALPHABET_SIZE;  i ++)
          if  (i != j)
              Ct [A + i] += 1.0;
#endif

        Div_By = 16;
        Mod_By = ALPHABET_SIZE;
        while  (Mod_By < MARKOV_LEN)
          {
           j = (Sub % Div_By) / Mod_By;
           A = Sub - j * Mod_By;
#if  PURINE_EDIT_DIST
           Ct [A + ((j + 2) % ALPHABET_SIZE) * Mod_By] += 1.0;
#else
           for  (i = 0;  i < ALPHABET_SIZE;  i ++)
             if  (i != j)
                 Ct [A + Mod_By * i] += 1.0;
#endif
           Div_By *= ALPHABET_SIZE;
           Mod_By *= ALPHABET_SIZE;
          }
       }

   return;
  }



int  Value  (char Ch)

/* Returns the "numeric" value of character  Ch . */

  {
   switch  (tolower (Ch))
     {
      case  'a' :
        return  0;
      case  'c' :
        return  1;
      case  'g' :
        return  2;
      case  't' :
        return  3;
      default :
//        fprintf (stderr, "ERROR:  Unexpected character '%c', assuming 'a'\n", Ch);
        return  0;
     }
  }
