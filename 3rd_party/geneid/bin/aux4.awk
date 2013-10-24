#USAGE= gawk -f aux.awk Background_P_file, Observed_P_file
BEGIN{
  while(getline<ARGV[1]>0){
      # read background probabilities
      BP[$1 $2]=$4;
      # Get conserved positions
      if ($4==1){
        ConservedPos[$2]=1;
        ConservedLet[$1 $2]=1;
      }
  }
  ARGV[1]="";
}

{
  if (ConservedPos[$2]==1){
    if (ConservedLet[$1 $2]==1)
      print $1, $2, 0;
    else
      print $1, $2, -9999;
  } else if ($4!=0 && BP[$1 $2]!=0 )
    print $1, $2, log($4/BP[$1 $2]);
  else print $1, $2, -9999;
}
