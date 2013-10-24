/************************************************************/
/* info.h                                                   */
/*                                                          */
/* Functions for calculating entropy, relative entropy,     */
/* conditional entropy, mutual information, conditional     */
/* mutual information, etc.                                 */
/*                                                          */
/************************************************************/

double I_mutual();
double entropy();
double rel_entropy();

/************************************************************/

/* EH[n] = Expected value of the sample entropy for a sample of n sequences	*/

double EH[26] = { 0.00000, 0.00000, 0.75000, 1.11090, 1.32399, 1.46291,
                           1.55923, 1.62900, 1.68129, 1.72155, 1.75328,
			   1.77879, 1.79966, 1.81699, 1.83159, 1.84403,
			   1.85475, 1.86408, 1.87227, 1.87952, 1.88598,
			   1.89177, 1.89699, 1.90172, 1.90604, 1.90998 };


/************************************************************/
 
double I_mutual(pab,pa,pb,na,nb) 
double *pab,*pa,*pb;
int na,nb;  
{ 
int i,j,ab; 
double papb,Iab=0.0;
 
for (i=0;i<na;++i) { 
  for (j=0;j<nb;++j) { ab=(na*i)+j;
    papb = pa[i]*pb[j];
    if (papb>0.0 && pab[ab]>0.0) Iab += pab[ab] * log(pab[ab]/papb)/log(2.0); 
    }  
  } 
 
return(Iab); 
 
}
 
/************************************************************/

double rel_entropy(a,b,n)	/* relative entropy of a rel to b */
double *a,*b;
int n;
{
int i;

double D=0.0;

for (i=0;i<n;++i)
  {
  if (a[i] > 0.0) D += a[i]*log(a[i]/b[i])/log(2.0);
  }

return(D);
}

/************************************************************/

double entropy(v,n)
double *v;
int n;
{
int i;
double H = 0.0;

for (i=0;i<n;++i) {
  if (v[i]<=0.0) continue;
  H += -1.0*v[i]*log(v[i])/log(2.0);
  }

return(H);
}

/************************************************************/
