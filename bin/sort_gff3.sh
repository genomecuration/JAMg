#!/bin/bash

awk 'BEGIN{OFS=FS="\t"} {
    if (NF == 0) next;
    if ($3=="gene") $3="1";
    else if ($3=="mRNA") $3="2";
    else if ($3=="exon") $3="3";
    else if ($3=="CDS") $3="4";
    else if ($3=="three_prime_UTR") $3="5";
    else if ($3=="five_prime_UTR") $3="6";
    print $0
}' $1 | sort --parallel=4 -k1,1 -k4,4n -k3,3n | awk 'BEGIN{OFS=FS="\t"} {if($3=="1") {print ""; $3="gene";} else if($3=="2") $3="mRNA"; else if($3=="3") $3="exon"; else if($3=="4") $3="CDS"; else if($3=="5") $3="three_prime_UTR"; else if($3=="6") $3="five_prime_UTR"; print $0}' > $1.sorted
