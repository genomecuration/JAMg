#! /usr/bin/gawk -f

BEGIN {
    begin=ARGV[1];
    end=ARGV[2];
    ARGV[1]=ARGV[2]="";
}
($2 >= begin && $2 <= end ){
    $2=$2-begin+1;print $2,$1,$3;
}
