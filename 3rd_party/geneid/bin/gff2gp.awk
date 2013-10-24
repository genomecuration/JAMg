#!/bin/gawk -f

# changes gff into golden path compliant format
# By changing the value of "context=" the user can select the number of bases of flanking sequence to be used.


{
  $4--;
  context=500; # can be changed according to how much sequence we want on 3' and 5' of gene
  gene=$9;
  seqID[gene]=$1;
  strand[gene]=$7;
  num_exons[gene]++;
  if (num_exons[gene]==1) {
    if ($4<1000) 
    {start_gene[gene]=1}
    else {
      start_gene[gene]=$4-context;
       } #else
    start_cds[gene]=$4;
    if (context<$4) {
    exons_starts[gene]=$4-context",";
    }
    else {    
      exons_starts[gene]=1",";
    }
    
} #if
  else {
    exons_starts[gene]=exons_starts[gene] $4",";
    exons_ends[gene]=exons_ends[gene] previous_end_gene[gene]"," 
      } #else
  previous_end_gene[gene]=$5;
  end_gene[gene]=$5+context;
  end_cds[gene]=$5;
  
}

END{
  for (gene in seqID)
    print gene, seqID[gene], strand[gene], start_gene[gene],end_gene[gene],start_cds[gene], end_cds[gene], num_exons[gene],exons_starts[gene],exons_ends[gene] end_gene[gene]"," ;


}
