//Copyright (c) 2003  by  Mihaela Pertea.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main  (int argc, char * argv [])
{

  double x;
  double y;

  x= strtod (argv[1],NULL);
  y=erf(x);

  printf("%g\n",y);
}
