#! /usr/bin/gawk -f
# frequency compute the frequency for every position of the sequences.
# gps, imim, sept 00
# USAGE: frequency k file
# K is the tuples number and the file must contain al the sequences with the same lenght

BEGIN {

  k=ARGV[1];
  ARGV[1]="";
  getline
      lenseq=length($2);

  alpha[1]="A";
  alpha[2]="C";
  alpha[3]="G";
  alpha[4]="T";

  PCOUNT=0;

  for (i=1;i<=4;i++)
     nx(1,k,alpha[i],MAT,PCOUNT,lenseq);
}

{
    for (i=1;i<=lenseq+1-k;i++) {
	if (substr($2,i,k) !~ /[^AGTC]/) {
	    MAT[substr($2,i,k),i]++;
	    MAT2[i]++;
	}
    }
}

END{
    for (i in MAT) {
	split (i,h,SUBSEP);
	printf "%s %2d %4d %.4f \n", h[1],h[2],MAT[i],MAT[i]/MAT2[h[2]] | "sort -k2,2n ";
	#printf "%s %2d %4d \n", h[1],h[2],MAT[i] | "sort +1n ";
    }
}

function nx(l,len,s,Mat,p,lenS,   i) {      # Funcio per inicialitzar les matrius
  if (l==len) {
      for (i=1;i<=(lenS+1-len);i++)
	  Mat[s,i]=p;
  }
  else {
     l++;
     for (i=1;i<=4;i++)
       nx(l,len,s alpha[i],Mat,p,lenS);
   }
}






