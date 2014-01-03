// Copyright (c) 1999, 2000, and 2001  by  Mihaela Pertea and A. L. Delcher.


//   A. L. Delcher
//
//     File:  ~delcher/TIGR/gene.h
//  Version:  1.01  31 Jul 97
//
//    Copyright (c) 1997 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  Common routines for DNA sequences.
//


#ifndef  __GENE_H_INCLUDED
#define  __GENE_H_INCLUDED

#define FALSE 0
#define TRUE 1

const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const int  MAX_LINE = 300;


char  Complement  (char);
char  Filter  (char);
int  Read_String  (FILE *, char * &, long int &, char [], int);
int  Read_Multi_String  (FILE *, char * &, long int &, char [], int,int);




#endif
