# MarkovMatrices.awk
# rgs, imim, april 98
# add pseudocounts, june 99
# usage=gawk -f MarkovMatrices order output_prefix cds_seqfile (.tbl)
#
# based on Borodosky and McIninch, Computers & Chemistry, 1993.

BEGIN {
  PCOUNT=5;

  k=ARGV[1];          #order of the markov chain
  pmat=ARGV[2];       # output matrix file name
  ARGV[1]=ARGV[2]="";

  # set pseudocounts 

  alpha[1]="A";
  alpha[2]="C";
  alpha[3]="G";
  alpha[4]="T";

  # for all k-tuples
  for (i=1;i<=4;i++) 
     nx(1,k,alpha[i],N0,4*PCOUNT);

  sizek=4^k;
  
  # for all k+1-tuples
  for (i=1;i<=4;i++) 
     nx(1,k+1,alpha[i],N,PCOUNT);

  sizek1=4^(k+1);

}
{
  sequence=toupper($2);
  if (sequence !~ /[^ACGT]/) { # consider only standard acgt

    lseq=length(sequence); 
    L+=(lseq-k);

    for (i=1;i<=lseq-k;i++) {
      ktuple=substr(sequence,i,k);
      ktuple1=substr(sequence,i,k+1);
      N0[ktuple]++; # intial probabilities
      N[ktuple1]++; # transition probabilities
    }
  }
}
END {
  # get number of k-tuples observed in each frame
  L0=L+sizek;

  for (t in N0) {       # initial probabilities
    P0[t]=N0[t]/(L0);
  }

  for (t in N) {        # transition probabilities
    tk=substr(t,1,k);   # k-tuple
    P[t,tk]=N[t]/N0[tk];
  }
    
  for (i in P0)
    print i,P0[i]  > pmat "." k ".initial";

  for (i in P) {
    split (i,x,SUBSEP);
    print x[1],x[2], P[i] > pmat "." k ".transition";
  }

}

function nx(l,len,s,Mat,p,   i) { 
  if (l==len) {
    Mat[s]=p;
    size++;
  }
  else {
     l++;
     for (i=1;i<=4;i++) 
       nx(l,len,s alpha[i], Mat,p);
   }
}
  
    
