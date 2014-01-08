use strict;
use GTF;

package Eval;

######################################################################
#
#  Definition of statistics and their data structures
#
######################################################################

# top level evaluate data structures
sub get_level_list{
    return qw (Gene
	       Transcript
	       Exon
               Nuc
	       Signal);
}

sub get_list_struct{
    my @level_list = get_level_list();
    my @gene_list = get_gene_type_list();
    my @tx_list = get_tx_type_list();
    my @exon_list = get_exon_type_list();
    my @nuc_list = get_nuc_type_list();
    my @signal_list = get_signal_type_list();
    my @gene_sub_list = get_gene_stat_list();
    my @tx_sub_list = get_tx_stat_list();
    my @exon_sub_list = get_exon_stat_list();
    my @nuc_sub_list = get_nuc_stat_list();
    my @signal_sub_list = get_signal_stat_list();
    my %struct = (Levels => \@level_list,
		  Gene => {Type => \@gene_list,
			   Stat => \@gene_sub_list},
		  Transcript => {Type => \@tx_list,
				 Stat => \@tx_sub_list},
		  Exon => {Type => \@exon_list,
			   Stat => \@exon_sub_list},
		  Nuc => {Type => \@nuc_list,
			  Stat => \@nuc_sub_list},
		  Signal => {Type => \@signal_list,
			     Stat => \@signal_sub_list});
    return %struct;
}

sub get_general_list_struct{
    my @level_list = get_level_list();
    my @gene_list = get_gene_type_list();
    my @tx_list = get_tx_type_list();
    my @exon_list = get_exon_type_list();
    my @nuc_list = get_nuc_type_list();
    my @signal_list = get_signal_type_list();
    my @gene_sub_list = get_gene_general_stat_list();
    my @tx_sub_list = get_tx_general_stat_list();
    my @exon_sub_list = get_exon_general_stat_list();
    my @nuc_sub_list = get_nuc_general_stat_list();
    my @signal_sub_list = get_signal_general_stat_list();
    my %struct = (Levels => \@level_list,
		  Gene => {Type => \@gene_list,
			   Stat => \@gene_sub_list},
		  Transcript => {Type => \@tx_list,
				 Stat => \@tx_sub_list},
		  Exon => {Type => \@exon_list,
			   Stat => \@exon_sub_list},
		  Nuc => {Type => \@nuc_list,
			  Stat => \@nuc_sub_list},
		  Signal => {Type => \@signal_list,
			     Stat => \@signal_sub_list});
    return %struct;
}

sub get_stats_struct{
    my %gene = _get_gene_type_struct();
    my %tx = _get_tx_type_struct();
    my %exon = _get_exon_type_struct();
    my %nuc = _get_nuc_type_struct();
    my %signal = _get_signal_type_struct();
    return (Gene => \%gene,
	    Exon => \%exon,
	    Transcript => \%tx,
	    Nuc => \%nuc,
	    Signal => \%signal);
}

# gene level evaluate data structures
sub get_gene_type_list{
    return qw (All);
}

sub get_gene_stat_list{
    my @list;
    my @genstat = get_gene_general_stat_list();
    foreach my $stat (@genstat){
	push @list, $stat;
    }
    my @add = get_gene_stat_type_list();
    my @substat = get_gene_substat_list();
    foreach my $add (@add){
	foreach my $substat (@substat){
	    push @list, $add."_".$substat;
	}
    }
    return @list;
}

sub get_gene_general_stat_list{
    return qw(Count
	      Ann_Count
	      Total_Transcripts
	      Transcripts_Per);
}

sub get_gene_stat_type_list{
    return qw(Consistent
	      Exact
	      Full_Exact
	      Genomic_Overlap
	      CDS_Overlap
	      All_Introns
	      Exact_UTR5
	      Exact_UTR3
	      Exact_Intron
	      Exact_Exon
	      Start_Codon
	      Stop_Codon
	      Start_Stop);
}

sub get_gene_substat_list{
    return qw(Pred_Count
	      Ann_Count
	      Specificity
	      Sensitivity);
}

sub _get_gene_type_hash{
    my @list = get_gene_type_list();
    my %gene;
    foreach my $type (@list){
	$gene{$type} = 0;
    }
    return %gene;
}

sub _get_gene_stat_hash{
    my @list = get_gene_stat_list();
    my %gene;
    foreach my $type (@list){
	$gene{$type} = 0;
    }
    return %gene;
}

sub _get_gene_type_struct{
    my @list = get_gene_type_list();
    my %gene;
    foreach my $type (@list){
	my %substruct = _get_gene_stat_hash(); 
	$gene{$type} = \%substruct;
    }
    return %gene;
}

# tx level evaluate data structures
sub get_tx_type_list{
    return qw (All 
	       Complete 
	       Incomplete
	       Single_Exon);
}

sub get_tx_stat_list{
    my @list;
    my @genstat = get_tx_general_stat_list();
    foreach my $stat (@genstat){
	push @list, $stat;
    }
    my @add = get_tx_stat_type_list();
    my @substat = get_tx_substat_list();
    foreach my $add (@add){
	foreach my $substat (@substat){
	    push @list, $add."_".$substat;
	}
    }
    return @list;
}

sub get_tx_general_stat_list{
    return qw(Count
	      Ann_Count
	      Average_Length
	      Median_Length
	      Total_Length
	      Average_Coding_Length
	      Median_Coding_Length
	      Total_Coding_Length
	      Average_Score
	      Total_Score
	      Ave_Exons_Per
		  Med_Exons_Per
	      Total_Exons);
}

sub get_tx_stat_type_list{
    return get_gene_stat_type_list();
}

sub get_tx_substat_list{
    return qw(Pred_Count
	      Ann_Count
	      Specificity
	      Sensitivity);
}

sub _get_tx_type_hash{
    my @list = get_tx_type_list();
    my %tx;
    foreach my $type (@list){
	$tx{$type} = 0;
    }
    return %tx;
}

sub _get_tx_stat_hash{
    my @list = get_tx_stat_list();
    my %tx;
    foreach my $type (@list){
	$tx{$type} = 0;
    }
    return %tx;
}

sub _get_tx_type_struct{
    my @list = get_tx_type_list();
    my %tx;
    foreach my $type (@list){
	my %substruct = _get_tx_stat_hash(); 
	$tx{$type} = \%substruct;
    }
    return %tx;
}

# exon level evaluate data structures
sub get_exon_type_list{
    return  qw(All
	       Initial
	       Internal
	       Terminal
	       Single
	       UTR3
	       UTR5
	       Intron
           InframeOptional
           FrameshiftOptional);
}

sub get_exon_stat_list{
    my @list;
    my @genstat = get_exon_general_stat_list();
    foreach my $stat (@genstat){
	push @list, $stat;
    }
    my @add = get_exon_stat_type_list();
    my @substat = get_exon_substat_list();
    foreach my $add (@add){
	foreach my $substat (@substat){
	    push @list, $add."_".$substat;
	}
    }
    return @list;
}

sub get_exon_general_stat_list{
    return qw(Count
	      Ann_Count
	      Average_Length
	      Median_Length
	      Total_Length
	      Average_Score
	      Total_Score);
}

sub get_exon_stat_type_list{
    return qw(Correct
	      Overlap
	      Overlap_80p
	      Splice_5
	      Splice_3);
}

sub get_exon_substat_list{
    return qw(Pred_Count
	      Ann_Count
	      Specificity
	      Sensitivity);
}

sub _get_exon_type_hash{
    my @list = get_exon_type_list();
    my %exon;
    foreach my $type (@list){
	$exon{$type} = 0;
    }
    return %exon;
}

sub _get_exon_stat_hash{
    my @list = get_exon_stat_list();
    my %exon;
    foreach my $type (@list){
	$exon{$type} = 0;
    }
    return %exon;
}

sub _get_exon_type_struct{
    my @list = get_exon_type_list();
    my %exon;
    foreach my $type (@list){	
	my %substruct = _get_exon_stat_hash(); 
	$exon{$type} = \%substruct;
    }
    return %exon;
}

# nuc level evaluate data structures
sub get_nuc_type_list{
    return get_exon_type_list();
}

sub get_nuc_stat_list{
    my @list;
    my @genstat = get_nuc_general_stat_list();
    foreach my $stat (@genstat){
	push @list, $stat;
    }
    my @add = get_nuc_stat_type_list();
    my @substat = get_nuc_substat_list();
    foreach my $add (@add){
	foreach my $substat (@substat){
	    push @list, $add."_".$substat;
	}
    }
    return @list;
}

sub get_nuc_general_stat_list{
    return qw(Count
	      Ann_Count);
}

sub get_nuc_stat_type_list{
    return qw (Correct);
}

sub get_nuc_substat_list{
    return qw(Pred_Count
	      Ann_Count
	      Specificity
	      Sensitivity);
}

sub _get_nuc_type_hash{
    my @list = get_nuc_type_list();
    my %nuc;
    foreach my $type (@list){
	$nuc{$type} = 0;
    }
    return %nuc;
}

sub _get_nuc_stat_hash{
    my @list = get_nuc_stat_list();
    my %nuc;
    foreach my $type (@list){
	$nuc{$type} = 0;
    }
    return %nuc;
}

sub _get_nuc_type_struct{
    my @list = get_nuc_type_list();
    my %nuc;
    foreach my $type (@list){	
	my %substruct = _get_nuc_stat_hash(); 
	$nuc{$type} = \%substruct;
    }
    return %nuc;
}

# signal level evaluate data structures
sub get_signal_type_list{
    return qw (Splice_Acceptor
	       Splice_Donor
	       Start_Codon
	       Stop_Codon);
}

sub get_signal_stat_list{
    my @list;
    my @genstat = get_signal_general_stat_list();
    foreach my $stat (@genstat){
	push @list, $stat;
    }
    my @add = get_signal_stat_type_list();
    my @substat = get_signal_substat_list();
    foreach my $add (@add){
	foreach my $substat (@substat){
	    push @list, $add."_".$substat;
	}
    }
    return @list;
}

sub get_signal_general_stat_list{
    return qw(Count
	      Ann_Count);
}

sub get_signal_stat_type_list{
    return qw (Correct);
}

sub get_signal_substat_list{
    return qw(Pred_Count
	      Ann_Count
	      Specificity
	      Sensitivity);
}

sub _get_signal_type_hash{
    my @list = get_signal_type_list();
    my %signal;
    foreach my $type (@list){
	$signal{$type} = 0;
    }    
    return %signal;
}

sub _get_signal_stat_hash{
    my @list = get_signal_stat_list();
    my %signal;
    foreach my $type (@list){
	$signal{$type} = 0;
    }
    return %signal;
}

sub _get_signal_type_struct{
    my @list = get_signal_type_list();
    my %signal;
    foreach my $type (@list){
	my %substruct = _get_signal_stat_hash(); 
	$signal{$type} = \%substruct;
    }
    return %signal;
}



######################################################################
#
#  evalutate($ann,$preds,$ann_name,$pred_names,$verbose)
#    Compares one set of annotation gtfs to one or more sets of 
#    prediction gtfs.
#
#    Parameters:
#      $ann     - reference to set(array) of annotation gtfs
#      $preds   - reference to array of sets(arrays) of prediction gtfs
#      $verbose - boolean turns verbose mode on/off
#
######################################################################
sub evaluate{
    my ($ann,$preds,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Evaluation Calculations\n";
    }
    my $start_time = time;
    my @a_genes;
    my @p_genes;
    my $a_data;
    my @p_data;
    foreach my $ann_gtf (@$ann){	
	my $genes = $ann_gtf->genes;
	push @a_genes, $genes;
    }
    $a_data = _get_stats(\@a_genes);
    _calculate_stats($a_data);
    foreach my $pred (@$preds){
	my $pred_genes;
	foreach my $pred_gtf (@$pred){
	    my $genes = $pred_gtf->genes;
	    push @$pred_genes, $genes;
	}
	push @p_genes, $pred_genes;
	my %stats_struct = get_stats_struct();
	push @p_data, \%stats_struct;
    }
    for(my $p = 0;$p <= $#p_genes;$p++){
	for(my $i = 0;$i <= $#a_genes;$i++){
	    compare_gene_lists($a_genes[$i],$p_genes[$p][$i],$p_data[$p]);
	}
	_calculate_stats($p_data[$p]);
    }
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    if($verbose){    
	print_time($total_time);
    }
    return ($a_data,@p_data);
}

# functions to initialize the data
sub _get_gene_tag{
    my %tag;
    my @types = get_gene_type_list();
    foreach my $type (@types){
	$tag{$type} = 0;
    }
    my @stats = get_gene_stat_type_list();
    foreach my $stat (@stats){
	$tag{$stat} = 0;
    }
    return \%tag;
}

sub _get_tx_tag{
    my %tag;
    my @types = get_tx_type_list();
    foreach my $type (@types){
	$tag{$type} = 0;
    }
    my @stats = get_tx_stat_type_list();
    foreach my $stat (@stats){
	$tag{$stat} = 0;
    }
    return \%tag;
}

sub _get_exon_tag{
    my %tag;
    my @types = get_exon_type_list();
    foreach my $type (@types){
	$tag{$type} = 0;
    }
    my @stats = get_exon_stat_type_list();
    foreach my $stat (@stats){
	$tag{$stat} = 0;
    }
    return \%tag;
}

sub _init_gene_tag{
    my ($gene) = @_;
    my $gene_tag = _get_gene_tag();
    $$gene_tag{All} = 1;
    $gene->set_tag($gene_tag);
    foreach my $tx (@{$gene->transcripts}){
	_init_tx_tag($tx);
    }
}

sub _init_tx_tag{
    my ($tx) = @_;
    my $tx_tag = _get_tx_tag();
    $$tx_tag{All} = 1;
    my $starts = $tx->start_codons;
    my $stops = $tx->stop_codons;
    my $cds = $tx->cds;
    if(($#$starts >= 0) && ($#$stops >= 0)){
	$$tx_tag{Complete} = 1;
	if($#$cds == 0){
	    $$tx_tag{Single_Exon} = 1;
	}
	else{
	    $$tx_tag{Single_Exon} = 0;
	}
    }
    else{
	$$tx_tag{Incomplete} = 1;	    
	$$tx_tag{Single_Exon} = 0;
    }
    $tx->set_tag($tx_tag);
    foreach my $cds (@{$tx->cds}){
	_init_exon_tag($cds);
    }
    foreach my $utr (@{$tx->utr3}){
	_init_utr3_tag($utr);
    }
    foreach my $utr (@{$tx->utr5}){
	_init_utr5_tag($utr);
    }
    foreach my $intron (@{$tx->introns}){
	_init_intron_tag($intron);
    }
}

sub _init_exon_tag{
    my ($cds) = @_;
    my $cds_tag = _get_exon_tag;
    $$cds_tag{All} = 1;
    my $type = $cds->subtype;
    $$cds_tag{$type} = 1;
    my $ase = $cds->ase;            # rpz add a tag for the alternative splicing
    $$cds_tag{$ase} = 1 if($ase);   # event if one exists.  Exons are put in
    $cds->set_tag($cds_tag);        # both buckets (subtype and ase).
}

sub _init_utr5_tag{
    my ($utr) = @_;
    my $u_tag = _get_exon_tag;
    $$u_tag{UTR5} = 1;
    $utr->set_tag($u_tag);
}

sub _init_utr3_tag{
    my ($utr) = @_;
    my $u_tag = _get_exon_tag;
    $$u_tag{UTR3} = 1;
    $utr->set_tag($u_tag);
}

sub _init_intron_tag{
    my ($intron) = @_;
    my $i_tag = _get_exon_tag;
    $$i_tag{Intron} = 1;
    $intron->set_tag($i_tag);
}

# free the memory
sub _clear_gene_tag{
    my ($gene) = @_;
    $gene->set_tag();
    foreach my $tx (@{$gene->transcripts}){
	_clear_tx_tag($tx);
    }
}

sub _clear_tx_tag{
    my ($tx) = @_;
    $tx->set_tag();
    foreach my $cds (@{$tx->cds}){
	_clear_exon_tag($cds);
    }
    foreach my $utr (@{$tx->utr5}){
	_clear_exon_tag($utr);
    }
    foreach my $utr (@{$tx->utr3}){
	_clear_exon_tag($utr);
    }
    foreach my $intron (@{$tx->introns}){
	_clear_exon_tag($intron);
    }
}

sub _clear_exon_tag{
    my ($exon) = @_;
    $exon->set_tag();
}

# functions to collect the data from the tags
sub _collect_gene_stats{    
    my ($genes,$data) = @_;
    my @all_txs;
    foreach my $gene (@$genes){
	my $info = $gene->tag;
	my $txs = $gene->transcripts;
	my @types = get_gene_type_list;
	my @stats = get_gene_stat_type_list;     
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Gene}{$type}{Count}++;
		$$data{Gene}{$type}{Total_Transcripts} += $#$txs + 1;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Gene}{$type}{$stat."_Pred_Count"}++;
		    }
		}
	    }
	}
	push @all_txs, @$txs;
    }
    @all_txs = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_txs;
    _collect_tx_stats(\@all_txs,$data);
}

sub _collect_ann_gene_stats{
    my ($genes,$data) = @_;
    my @all_txs;
    foreach my $gene (@$genes){
	my $info = $gene->tag;
	my $txs = $gene->transcripts;
	my @types = get_gene_type_list;
	my @stats = get_gene_stat_type_list;     
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Gene}{$type}{Ann_Count}++;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Gene}{$type}{$stat."_Ann_Count"}++;
		    }
		}
	    }
	}
	push @all_txs, @$txs;
    }
    @all_txs = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_txs;
    _collect_ann_tx_stats(\@all_txs,$data);
}

sub _collect_tx_stats{
    my ($txs,$data) = @_;
    my @all_cds;
    my @all_introns;
    my @all_utr5;
    my @all_utr3;
    foreach my $tx (@$txs){
	my $info = $tx->tag;
	my $cds = $tx->cds;
	my $introns = $tx->introns;
	my $len = $tx->length;
	my $clen = $tx->coding_length;
	my $score = $tx->score;
	my @types = get_tx_type_list;
	my @stats = get_tx_stat_type_list;     
	my $starts = $tx->start_codons;
	my $stops = $tx->stop_codons;
	if($#$starts >= 0){
	    $$data{Signal}{Start_Codon}{Count}++;
	    if($$info{Start_Codon}){
		$$data{Signal}{Start_Codon}{Correct_Pred_Count}++;
	    }
	}
	if($#$stops >= 0){
	    $$data{Signal}{Stop_Codon}{Count}++;
	    if($$info{Stop_Codon}){
		$$data{Signal}{Stop_Codon}{Correct_Pred_Count}++;
	    }
	}
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Transcript}{$type}{Count}++;
		$$data{Transcript}{$type}{Total_Length} += $len;
		push @{$$data{Transcript}{$type}{Length_Array}}, $len;
		$$data{Transcript}{$type}{Total_Coding_Length} += $clen;
		push @{$$data{Transcript}{$type}{Coding_Length_Array}}, $clen;
		$$data{Transcript}{$type}{Total_Exons} += $#$cds + 1;
		push @{$$data{Transcript}{$type}{Exons_Per_Array}}, $#$cds + 1;
		$$data{Transcript}{$type}{Total_Score} += $score;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Transcript}{$type}{$stat."_Pred_Count"}++;
		    }
		}
	    }
	}
	push @all_cds, @$cds;
	push @all_introns, @$introns;
	push @all_utr5, @{$tx->utr5};
	push @all_utr3, @{$tx->utr3};
    }
    
    @all_cds = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_cds;
    @all_utr5 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr5;
    @all_utr3 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr3;
    @all_introns = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_introns;

    _collect_exon_stats(\@all_cds,$data);
    _collect_exon_stats(\@all_utr5,$data);
    _collect_exon_stats(\@all_utr3,$data);
    _collect_exon_stats(\@all_introns,$data);
}

sub _collect_ann_tx_stats{
    my ($txs,$data) = @_;
    my @all_cds;
    my @all_introns;
    my @all_utr5;
    my @all_utr3;
    foreach my $tx (@$txs){    
	my $info = $tx->tag;
	my $cds = $tx->cds;
	my $introns = $tx->introns;
	my $len = $tx->length;
	my $clen = $tx->coding_length;
	my $score = $tx->score;
	my @types = get_tx_type_list;
	my @stats = get_tx_stat_type_list;     
	my $starts = $tx->start_codons;
	my $stops = $tx->stop_codons;
	if($#$starts >= 0){
	    $$data{Signal}{Start_Codon}{Ann_Count}++;
	    if($$info{Start_Codon}){
		$$data{Signal}{Start_Codon}{Correct_Ann_Count}++;
	    }
	}
	if($#$stops >= 0){
	    $$data{Signal}{Stop_Codon}{Ann_Count}++;
	    if($$info{Stop_Codon}){
		$$data{Signal}{Stop_Codon}{Correct_Ann_Count}++;
	    }
	}
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Transcript}{$type}{Ann_Count}++;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Transcript}{$type}{$stat."_Ann_Count"}++;
		    }
		}
	    }
	}
	push @all_cds, @$cds;
	push @all_utr5, @{$tx->utr5};
	push @all_utr3, @{$tx->utr3};
	push @all_introns, @$introns;
    }
    
    @all_cds = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_cds;
    @all_utr5 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr5;
    @all_utr3 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr3;
    @all_introns = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_introns;
    _collect_ann_exon_stats(\@all_cds,$data);
    _collect_ann_exon_stats(\@all_utr5,$data);
    _collect_ann_exon_stats(\@all_utr3,$data);
    _collect_ann_exon_stats(\@all_introns,$data);
}

sub _collect_exon_stats{
    my ($exons,$data) = @_; 
    my @types = get_exon_type_list;
    my @stats = get_exon_stat_type_list;         
    my @skip;
    my %all_nucs;
    my %match_nucs;
    foreach my $type (@types){
	$all_nucs{$type} = [];
	$match_nucs{$type} = [];
    }
    for(my $i = 0;$i <= $#$exons;$i++){
	my $old_info = $$exons[$i]->tag;
	foreach my $type (@types){
	    if($$old_info{$type}){
		push @{$all_nucs{$type}}, [$$exons[$i]->start,$$exons[$i]->stop];
		if(defined($$old_info{Correct_Nucs})){
		    push @{$match_nucs{$type}}, @{$$old_info{Correct_Nucs}};	
		}
	    }
	}
	if($$exons[$i]->equals($$exons[$i+1])){
	    my $new_info = $$exons[$i+1]->tag;
	    foreach my $key (keys %$old_info){
		unless($key eq "Correct_Nucs"){
		    if($$old_info{$key}){
			$$new_info{$key} = 1;
		    }
		}
	    }
	    $skip[$i] = 1;
	}
	else{
	    $skip[$i] = 0;
	}
    }
    foreach my $type (@types){
	$match_nucs{$type} = [sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @{$match_nucs{$type}}];
	$all_nucs{$type} = [sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @{$all_nucs{$type}}];
	my $all_count = 0;
	my $match_count = 0;
	my $last_all = -1;
	my $last_match = -1;
	foreach my $n (@{$all_nucs{$type}}){
	    if($$n[1] > $last_all){
		$all_count += $$n[1] - max($$n[0],$last_all+1) + 1;		
		$last_all = $$n[1];
	    }
	}
	foreach my $n (@{$match_nucs{$type}}){
	    if($$n[1] > $last_match){
		$match_count += $$n[1] - max($$n[0],$last_match+1) + 1;
		$last_match = $$n[1];
	    }
	}
	$$data{Nuc}{$type}{Count} += $all_count;
	$$data{Nuc}{$type}{Correct_Pred_Count} += $match_count;
    }
    for(my $i = 0;$i <= $#$exons;$i++){
	my $exon = $$exons[$i];
	if($skip[$i]){
	    next;
	}
	my $info = $exon->tag;
	my $len = $exon->length;    
	my $score = $exon->score;    
	my $cnucs = 0;
	foreach my $interval (@{$$info{Correct_Nucs}}){
	    $cnucs += $$interval[1] - $$interval[0] + 1;
	}
	if($$info{Initial} || $$info{Internal}){
	    $$data{Signal}{Splice_Donor}{Count}++;
	    if($$info{Splice_3}){
		$$data{Signal}{Splice_Donor}{Correct_Pred_Count}++;
	    }
	}
	if($$info{Terminal} || $$info{Internal}){
	    $$data{Signal}{Splice_Acceptor}{Count}++;	    
	    if($$info{Splice_5}){
		$$data{Signal}{Splice_Acceptor}{Correct_Pred_Count}++;
	    }
	}
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Exon}{$type}{Count}++;
		$$data{Exon}{$type}{Total_Length} += $len;
		push @{$$data{Exon}{$type}{Length_Array}}, $len;
		$$data{Exon}{$type}{Total_Score} += $score;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Exon}{$type}{$stat."_Pred_Count"}++;
		    }
		}
	    }	
	}
    }
}

sub _collect_ann_exon_stats{
    my ($exons,$data) = @_;
    my @types = get_exon_type_list;
    my @stats = get_exon_stat_type_list;     
    my @skip;
    my %all_nucs;
    my %match_nucs;
    foreach my $type (@types){
	$all_nucs{$type} = [];
	$match_nucs{$type} = [];
    }
    for(my $i = 0;$i <= $#$exons;$i++){
	my $old_info = $$exons[$i]->tag;
	foreach my $type (@types){
	    if($$old_info{$type}){
		push @{$all_nucs{$type}}, [$$exons[$i]->start,$$exons[$i]->stop];
		if(defined($$old_info{Correct_Nucs})){
		    push @{$match_nucs{$type}}, @{$$old_info{Correct_Nucs}};	
		}
	    }
	}
	if($$exons[$i]->equals($$exons[$i+1])){
	    my $new_info = $$exons[$i+1]->tag;
	    foreach my $key (keys %$old_info){
		unless($key eq "Correct_Nucs"){
		    if($$old_info{$key}){
			$$new_info{$key} = 1;
		    }
		}
	    }
	    $skip[$i] = 1;
	}
	else{
	    $skip[$i] = 0;
	}
    }
    foreach my $type (@types){
	$match_nucs{$type} = [sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @{$match_nucs{$type}}];
	$all_nucs{$type} = [sort {$$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @{$all_nucs{$type}}];
	my $all_count = 0;
	my $match_count = 0;
	my $last_all = -1;
	my $last_match = -1;
	foreach my $n (@{$all_nucs{$type}}){
	    if($$n[1] > $last_all){
		$all_count += $$n[1] - max($$n[0],$last_all+1) + 1;		
		$last_all = $$n[1];
	    }
	}
	foreach my $n (@{$match_nucs{$type}}){
	    if($$n[1] > $last_match){
		$match_count += $$n[1] - max($$n[0],$last_match+1) + 1;
		$last_match = $$n[1];
	    }
	}
	$$data{Nuc}{$type}{Ann_Count} += $all_count;
	$$data{Nuc}{$type}{Correct_Ann_Count} += $match_count;
    }
    for(my $i = 0;$i <= $#$exons;$i++){
	my $exon = $$exons[$i];
	if($skip[$i]){
	    next;
	}
	my $info = $exon->tag;
	my $len = $exon->length;    
	my $score = $exon->score;    
	my $cnucs = 0;
	foreach my $interval (@{$$info{Correct_Nucs}}){
	    $cnucs += $$interval[1] - $$interval[0] + 1;
	}
	if($$info{Initial} || $$info{Internal}){
	    $$data{Signal}{Splice_Donor}{Ann_Count}++;
	    if($$info{Splice_3}){
		$$data{Signal}{Splice_Donor}{Correct_Ann_Count}++;
	    }
	}
	if($$info{Terminal} || $$info{Internal}){
	    $$data{Signal}{Splice_Acceptor}{Ann_Count}++;	    
	    if($$info{Splice_5}){
		$$data{Signal}{Splice_Acceptor}{Correct_Ann_Count}++;
	    }
	}
	foreach my $type (@types){
	    if($$info{$type}){
		$$data{Exon}{$type}{Ann_Count}++;
		foreach my $stat (@stats){
		    if($$info{$stat}){
			$$data{Exon}{$type}{$stat."_Ann_Count"}++;
		    }
		}
	    }
	}
    }
}

# function to compute the stats
sub _calculate_stats{
    my ($data) = @_;
    _calculate_gene_stats($data);
    _calculate_tx_stats($data);
    _calculate_exon_stats($data);
    _calculate_nuc_stats($data);
    _calculate_signal_stats($data);
}

sub _calculate_gene_stats{
    my ($data) = @_;
    my @types = get_gene_type_list;
    my @stats = get_gene_stat_type_list;
    foreach my $type (@types){
	if($$data{Gene}{$type}{Count} > 0){
	    $$data{Gene}{$type}{Transcripts_Per} =  
		$$data{Gene}{$type}{Total_Transcripts} / 
		    $$data{Gene}{$type}{Count};
	}
	foreach my $stat (@stats){
	    if($$data{Gene}{$type}{Count} > 0){
		$$data{Gene}{$type}{$stat."_Specificity"} =
		    100 * $$data{Gene}{$type}{$stat."_Pred_Count"} / 
			$$data{Gene}{$type}{Count};
	    }
	    if($$data{Gene}{$type}{Ann_Count} > 0){
		$$data{Gene}{$type}{$stat."_Sensitivity"} =
		    100 * $$data{Gene}{$type}{$stat."_Ann_Count"} / 
			$$data{Gene}{$type}{Ann_Count};
	    }   
	}
    }
}

sub _calculate_tx_stats{
    my ($data) = @_;
    my @types = get_tx_type_list;
    my @stats = get_tx_stat_type_list;
	my @sorted;
    foreach my $type (@types){
	if($$data{Transcript}{$type}{Count} > 0){
	    $$data{Transcript}{$type}{Average_Score} =  
		$$data{Transcript}{$type}{Total_Score} / 
		    $$data{Transcript}{$type}{Count};
	    $$data{Transcript}{$type}{Average_Length} =  
		$$data{Transcript}{$type}{Total_Length} / 
		    $$data{Transcript}{$type}{Count};
		@sorted = sort {$a <=> $b} @{$$data{Transcript}{$type}{Length_Array}}; 
	    $$data{Transcript}{$type}{Median_Length} =  
		    $sorted[$$data{Transcript}{$type}{Count}/2];
	    $$data{Transcript}{$type}{Average_Coding_Length} =  
		$$data{Transcript}{$type}{Total_Coding_Length} / 
		    $$data{Transcript}{$type}{Count};
		@sorted = sort {$a <=> $b} @{$$data{Transcript}{$type}{Coding_Length_Array}}; 
	    $$data{Transcript}{$type}{Median_Coding_Length} =  
		    $sorted[$$data{Transcript}{$type}{Count}/2];
	    $$data{Transcript}{$type}{Ave_Exons_Per} =
		$$data{Transcript}{$type}{Total_Exons} / 
		    $$data{Transcript}{$type}{Count};
		@sorted = sort {$a <=> $b} @{$$data{Transcript}{$type}{Exons_Per_Array}};
		$$data{Transcript}{$type}{Med_Exons_Per} =
			$sorted[$$data{Transcript}{$type}{Count}/2];
	}
	foreach my $stat (@stats){
	    if($$data{Transcript}{$type}{Count} > 0){
		$$data{Transcript}{$type}{$stat."_Specificity"} =
		    100 * $$data{Transcript}{$type}{$stat."_Pred_Count"} / 
			$$data{Transcript}{$type}{Count};
	    }
	    if($$data{Transcript}{$type}{Ann_Count} > 0){
		$$data{Transcript}{$type}{$stat."_Sensitivity"} =
		    100 * $$data{Transcript}{$type}{$stat."_Ann_Count"} / 
			$$data{Transcript}{$type}{Ann_Count};
	    }   
	}
    }
}

sub _calculate_exon_stats{
    my ($data) = @_;
    my @types = get_exon_type_list;
    my @stats = get_exon_stat_type_list;
	my @sorted;
    foreach my $type (@types){
	if($$data{Exon}{$type}{Count}){
	    $$data{Exon}{$type}{Average_Length} = 
		$$data{Exon}{$type}{Total_Length} / 
		    $$data{Exon}{$type}{Count};
		@sorted = sort {$a <=> $b} @{$$data{Exon}{$type}{Length_Array}}; 
	    $$data{Exon}{$type}{Median_Length} = 
		    $sorted[$$data{Exon}{$type}{Count}/2];
	    $$data{Exon}{$type}{Average_Score} = 
		$$data{Exon}{$type}{Total_Score} / 
		    $$data{Exon}{$type}{Count};	    
	}
	foreach my $stat (@stats){
	    if($$data{Exon}{$type}{Count}){
		$$data{Exon}{$type}{$stat."_Specificity"} =
		    100 * $$data{Exon}{$type}{$stat."_Pred_Count"} / 
			$$data{Exon}{$type}{Count};
	    }
	    if($$data{Exon}{$type}{Ann_Count}){
		$$data{Exon}{$type}{$stat."_Sensitivity"} =
		    100 * $$data{Exon}{$type}{$stat."_Ann_Count"} / 
			$$data{Exon}{$type}{Ann_Count};
	    }
	}
    }
}

sub _calculate_nuc_stats{
    my ($data) = @_;
    my @types = get_nuc_type_list;
    my @stats = get_nuc_stat_type_list;
    foreach my $type (@types){
	foreach my $stat (@stats){
	    if($$data{Nuc}{$type}{Count}){
		$$data{Nuc}{$type}{$stat."_Specificity"} =
		    100 * $$data{Nuc}{$type}{$stat."_Pred_Count"} / 
			$$data{Nuc}{$type}{Count};
	    }
	    if($$data{Nuc}{$type}{Ann_Count}){
		$$data{Nuc}{$type}{$stat."_Sensitivity"} =
		    100 * $$data{Nuc}{$type}{$stat."_Ann_Count"} / 
			$$data{Nuc}{$type}{Ann_Count};
	    }
	}
    }
}

sub _calculate_signal_stats{
    my ($data) = @_;
    my @types = get_signal_type_list;
    my @stats = get_signal_stat_type_list;
    foreach my $type (@types){
	foreach my $stat (@stats){
	    if($$data{Signal}{$type}{Count}){
		$$data{Signal}{$type}{$stat."_Specificity"} =
		    100 * $$data{Signal}{$type}{$stat."_Pred_Count"} / 
			$$data{Signal}{$type}{Count};
	    }
	    if($$data{Signal}{$type}{Ann_Count}){
		$$data{Signal}{$type}{$stat."_Sensitivity"} =
		    100 * $$data{Signal}{$type}{$stat."_Ann_Count"} / 
			$$data{Signal}{$type}{Ann_Count};
	    }
	}
    }
}

# functions to compare lists of objects
sub _compare_object_lists{
    my ($anns,$preds,$data,$init_func,$clear_func,$compare_func,
	$collect_func,$ann_collect_func) = @_;
    my $next_a = 0;
    my $a = 0;
    my $max_init = -1;
    my $last_cleared = -1;
    my @pending_preds;
    my @pending_anns;
    my $last_pred_stop = -1;
    my $last_ann_stop = -1;    
    if($#$anns >= 0){
	$last_ann_stop = $$anns[0]->stop;
    }
    if($#$preds >= 0){
	$last_pred_stop = $$preds[0]->stop;
    }
    for(my $p = 0;$p <= $#$preds;$p++){
	&$init_func($$preds[$p]);
	my $p_start = $$preds[$p]->start;
	my $p_stop = $$preds[$p]->stop;
	my $p_strand = $$preds[$p]->strand;	
	while(($next_a <= $#$anns) && ($p_start > $$anns[$next_a]->stop)){
	    if($next_a > $max_init){
		&$init_func($$anns[$next_a]);
		$max_init = $next_a;
	    }
	    $next_a++;	
	}
	my $a = $next_a;
	while(($a <= $#$anns) && ($p_stop >= $$anns[$a]->start)){
	    if($a > $max_init){
		&$init_func($$anns[$a]);
		$max_init = $a;
	    }
	    &$compare_func($$anns[$a],$$preds[$p]);
	    $a++;
	}
	if($$preds[$p]->start > $last_pred_stop){
	    &$collect_func(\@pending_preds,$data);
	    foreach my $pred (@pending_preds){
		&$clear_func($pred);
	    }
	    @pending_preds = ();
	    $last_pred_stop = 0;
	}
	if($last_pred_stop < $$preds[$p]->stop){
	    $last_pred_stop = $$preds[$p]->stop;
	}
	push @pending_preds, $$preds[$p];
	for(my $i = $last_cleared + 1;$i < $next_a;$i++){
	    if($$anns[$i]->start > $last_ann_stop){
		&$ann_collect_func(\@pending_anns,$data);
		foreach my $ann (@pending_anns){
		    &$clear_func($ann);
		}
		@pending_anns = ();
		$last_ann_stop = 0;
	    }
	    if($last_ann_stop < $$anns[$i]->stop){
		$last_ann_stop = $$anns[$i]->stop;
	    }
	    push @pending_anns, $$anns[$i];
	    $last_cleared = $i;
	}
    }   
    &$collect_func(\@pending_preds,$data);
    for(my $i = $max_init + 1;$i <= $#$anns;$i++){
	&$init_func($$anns[$i]);	
    }
    for(my $i = $last_cleared + 1;$i <= $#$anns;$i++){
	push @pending_anns, $$anns[$i];
    }
    &$ann_collect_func(\@pending_anns,$data);
    foreach my $ann (@pending_anns){
	&$clear_func($ann);
    }
}

sub compare_gene_lists{
    my ($anns,$preds,$data) = @_;
    _compare_object_lists($anns,$preds,$data,\&_init_gene_tag,\&_clear_gene_tag,
		          \&_compare_genes,\&_collect_gene_stats,
			  \&_collect_ann_gene_stats);
}

sub compare_tx_lists{
    my ($anns,$preds,$data) = @_;
    _compare_object_lists($anns,$preds,$data,\&_init_tx_tag,\&_clear_tx_tag,
			  \&_compare_txs,\&_collect_tx_stats,\&_collect_ann_tx_stats);
}

sub compare_exon_lists{
    my ($anns,$preds,$data) = @_;
    _compare_object_lists($anns,$preds,$data,\&_init_exon_tag,\&_clear_exon_tag,
			  \&_compare_features,\&_collect_exon_stats,
			  \&_collect_ann_exon_stats);
}

# functions to compare two objects
sub _compare_genes{
    my ($a_gene,$p_gene) = @_;
    if($p_gene->strand ne $a_gene->strand){
	return;
    }		
    my $a_txs = $a_gene->transcripts;
    my $p_txs = $p_gene->transcripts;
    my $next_a = 0;
    for(my $p = 0;$p <= $#$p_txs;$p++){
	my $p_start = $$p_txs[$p]->start;
	my $p_stop = $$p_txs[$p]->stop;
	while(($next_a <= $#$a_txs) && ($p_start > $$a_txs[$next_a]->stop)){
	    $next_a++;
	}
	my $a = $next_a;
	while(($a <= $#$a_txs) && ($p_stop >= $$a_txs[$a]->start)){
	    _compare_txs($$a_txs[$a],$$p_txs[$p]);
	    $a++;
	}
    }
}

sub _compare_txs{
    my ($a_tx,$p_tx) = @_;
    my $a_gene = $a_tx->gene;
    my $p_gene = $p_tx->gene;
    my $ag_info = $a_gene->tag;
    my $pg_info = $p_gene->tag;
    my $a_info = $a_tx->tag;
    my $p_info = $p_tx->tag;
    my $ann_starts = $a_tx->start_codons;
    my $ann_stops = $a_tx->stop_codons;
    my $ann_cds = $a_tx->cds;
    my $ann_introns = $a_tx->introns;
    my $ann_utr5 = $a_tx->utr5;
    my $ann_utr3 = $a_tx->utr3;
    my $pred_starts = $p_tx->start_codons;
    my $pred_stops = $p_tx->stop_codons;
    my $pred_cds = $p_tx->cds;
    my $pred_introns = $p_tx->introns;
    my $pred_utr5 = $p_tx->utr5;
    my $pred_utr3 = $p_tx->utr3;
    my $consistent = 1;
    my $exact = 1;
    my $start_codon = 0;
    my $stop_codon = 0;
    my $pred_complete = 1;
    my $ann_complete = 1;
    if($p_tx->strand ne $a_tx->strand){
	return;
    }
    #check that pred overlaps all of ann
    if($p_tx->coding_start > $a_tx->coding_start){
	$consistent = 0;
	$exact = 0;
    }
    if($p_tx->coding_stop < $a_tx->coding_stop){
	$consistent = 0;
	$exact = 0;
    }
    if(($a_tx->stop > $p_tx->start) &&
       ($a_tx->start < $p_tx->stop)){
	$$a_info{Genomic_Overlap} = 1;
	$$p_info{Genomic_Overlap} = 1;
	$$ag_info{Genomic_Overlap} = 1;
	$$pg_info{Genomic_Overlap} = 1;
    }
    else{
	return;
    }
    #check start and stop codons match
    if(($#$ann_starts >= 0) && ($#$pred_starts >= 0)){
	if($#$ann_starts != $#$pred_starts){
	    $consistent = 0;
	    $exact = 0;
	}
	else{
	    $start_codon = 1;
	    for(my $i = 0;$i <= $#$ann_starts;$i++){
		if(($$ann_starts[$i]->start != $$pred_starts[$i]->start) ||
		   ($$ann_starts[$i]->stop != $$pred_starts[$i]->stop)){
		    $start_codon = 0;
		    $consistent = 0;
		    $exact = 0;
		}
	    }
	}
    }
    elsif(($#$ann_starts >= 0) && ($#$pred_starts == -1)){
	$pred_complete = 0;
	$consistent = 0;
	$exact = 0;
    }
    elsif(($#$ann_starts == -1) && ($#$pred_starts >= 0)){
	$ann_complete = 0;
	$exact = 0;
    }
    else{
	$pred_complete = 0;
	$ann_complete = 0;
    }
    if(($#$ann_stops >= 0) && ($#$pred_stops >= 0)){
	if($#$ann_stops != $#$pred_stops){
	    $consistent = 0;
	    $exact = 0;
	}
	else{
	    $stop_codon = 1;
	    for(my $i = 0;$i <= $#$ann_stops;$i++){
		if(($$ann_stops[$i]->start != $$pred_stops[$i]->start) ||
		   ($$ann_stops[$i]->stop != $$pred_stops[$i]->stop)){
		    $stop_codon = 0;
		    $consistent = 0;
		    $exact = 0;
		}
	    }
	}
    }
    elsif(($#$ann_stops >= 0) && ($#$pred_stops == -1)){
	$pred_complete = 0;
	$consistent = 0;
	$exact = 0;
    }
    elsif(($#$ann_stops == -1) && ($#$pred_stops >= 0)){
	$ann_complete = 0;
	$exact = 0;
    }
    else{
	$pred_complete = 0;
	$ann_complete = 0;
    }
    if($start_codon && $stop_codon){
	$$p_info{Start_Stop} = 1;
	$$a_info{Start_Stop} = 1;
	$$pg_info{Start_Stop} = 1;
	$$ag_info{Start_Stop} = 1;
    }
    if($start_codon){
	$$p_info{Start_Codon} = 1;
	$$a_info{Start_Codon} = 1;
	$$pg_info{Start_Codon} = 1;
	$$ag_info{Start_Codon} = 1;
    }
    if($stop_codon){
	$$p_info{Stop_Codon} = 1;
	$$a_info{Stop_Codon} = 1;
	$$pg_info{Stop_Codon} = 1;
	$$ag_info{Stop_Codon} = 1;
    }
    
    my $ann_low = -1;
    my $ann_high = -1;
    if($a_tx->strand eq '+'){
	if($#$ann_starts >= 0){
	    $ann_low = $$ann_starts[0]->start;
	}
	if($#$ann_stops >= 0){
	    $ann_high = $$ann_stops[$#$ann_stops]->stop;
	}
    }
    else{
	if($#$ann_starts >= 0){
	    $ann_high = $$ann_starts[$#$ann_starts]->stop;
	}
	if($#$ann_stops >= 0){
	    $ann_low = $$ann_stops[0]->start;
	}
    }
    #check introns
    my $all_introns = 1;
    my $exact_intron = 0;
    if($#$ann_introns == -1){
	$all_introns = 0;
    }
    my $next_p = 0;
    if($#$ann_introns >= 0){
	while(($next_p <= $#$pred_introns) && 
	      ($$ann_introns[0]->start > $$pred_introns[$next_p]->stop)){
	    $next_p++;
	}
	if($next_p > 0){
	    if($ann_low >= 0){
		$consistent = 0;
		$all_introns = 0;
	    }
	    $exact = 0;
	}
    }
    else{
	$all_introns = 0;
    }
    my $max_checked = 0;
    for(my $a = 0;$a <= $#$ann_introns;$a++){
	my $a_start = $$ann_introns[$a]->start;
	my $a_stop = $$ann_introns[$a]->stop;
	while(($next_p <= $#$pred_introns) && 
	      ($a_start > $$pred_introns[$next_p]->stop)){
	    $next_p++;
	    if($next_p > ($max_checked + 1)){
		$all_introns = 0;
		$consistent = 0;
		$exact = 0;
	    }
	}
	my $p = $next_p;
	my $match = 0;
	# if its a match only one will match since features can't overlap
	while(($p <= $#$pred_introns) && ($a_stop >= $$pred_introns[$p]->start)){
	    $match = _compare_features($$ann_introns[$a],$$pred_introns[$p]);
	    $max_checked = $p;
	    $p++;
	}
	if($match){
	    $exact_intron = 1;
	}
	else{
	    $all_introns = 0;
	    $consistent = 0;
	    $exact = 0;
	}
    }
    if($max_checked < $#$pred_introns){
	$exact = 0;
	if($ann_high >= 0){
	    $consistent = 0;
	    $all_introns = 0;
	}
    }
    if($all_introns){
	$$a_info{All_Introns} = 1;
	$$p_info{All_Introns} = 1;
	$$ag_info{All_Introns} = 1;
	$$pg_info{All_Introns} = 1;
    }
    if($exact_intron){
    	$$a_info{Exact_Intron} = 1;
	$$p_info{Exact_Intron} = 1;
    	$$ag_info{Exact_Intron} = 1;
	$$pg_info{Exact_Intron} = 1;
    }

    #check utr5
    my $utr5_match = 0;
    my $a = 0;
    for(my $p = 0;$p <= $#$pred_utr5;$p++){
	my $p_start = $$pred_utr5[$p]->start;
	my $p_stop = $$pred_utr5[$p]->stop;
	my $p_strand = $$pred_utr5[$p]->strand;	
	my $old_a = $a;
	# features can't overlap others in the same tx and must 
	# be on the same strand with everything we dealing with here
	while(($a <= $#$ann_utr5) && ($p_start > $$ann_utr5[$a]->stop)){
	    $a++;
	}
	while(($a <= $#$ann_utr5) && ($p_stop >= $$ann_utr5[$a]->start)){
	    $utr5_match += _compare_features($$ann_utr5[$a],$$pred_utr5[$p]);
	    $a++;
	}
	# an ann can overlap multiple preds (but no anns) so this makes
	# sure that it gets checked against them all
	if($a > $old_a){
	    $a--;
	}
    }
    if($start_codon && $#$pred_utr5 == $#$ann_utr5 && ($#$pred_utr5 + 1 == $utr5_match)){
	$$a_info{Exact_UTR5} = 1;
	$$p_info{Exact_UTR5} = 1;
	$$ag_info{Exact_UTR5} = 1;
	$$pg_info{Exact_UTR5} = 1;
	$utr5_match = 1;
    }
    else{
	$utr5_match = 0;
    }
    #check utr3
    my $utr3_match = 0;
    $a = 0;
    for(my $p = 0;$p <= $#$pred_utr3;$p++){
	my $p_start = $$pred_utr3[$p]->start;
	my $p_stop = $$pred_utr3[$p]->stop;
	my $p_strand = $$pred_utr3[$p]->strand;	
	my $old_a = $a;
	# features can't overlap others in the same tx and must 
	# be on the same strand with everything we dealing with here
	while(($a <= $#$ann_utr3) && ($p_start > $$ann_utr3[$a]->stop)){
	    $a++;
	}
	while(($a <= $#$ann_utr3) && ($p_stop >= $$ann_utr3[$a]->start)){
	    $utr3_match += _compare_features($$ann_utr3[$a],$$pred_utr3[$p]);
	    $a++;
	}
	# an ann can overlap multiple preds (but no anns) so this makes
	# sure that it gets checked against them all
	if($a > $old_a){
	    $a--;
	}
    }
    if($stop_codon && $#$pred_utr3 == $#$ann_utr3 && ($#$pred_utr3 + 1 == $utr3_match)){
	$$a_info{Exact_UTR3} = 1;
	$$p_info{Exact_UTR3} = 1;
	$$ag_info{Exact_UTR3} = 1;
	$$pg_info{Exact_UTR3} = 1;
	$utr3_match = 1;
    }
    else{
	$utr3_match = 0;
    }

    #check cds
    my $exact_cds = 0;
    my $nuc_overlap = 0;
    $next_p = 0;
    if($#$ann_cds >= 0){
	while(($next_p <= $#$pred_cds) && 
	      ($$ann_cds[0]->start > $$pred_cds[$next_p]->stop)){
	    $next_p++;
	}
	if($next_p > 0){
	    if($ann_low >= 0){
		$consistent = 0;
	    }
	    $exact = 0;
	}
    }
    $max_checked = 0;
    for(my $a = 0;$a <= $#$ann_cds;$a++){
	my $a_start = $$ann_cds[$a]->start;
	my $a_stop = $$ann_cds[$a]->stop;
	while(($next_p <= $#$pred_cds) && 
	      ($a_start > $$pred_cds[$next_p]->stop)){
	    $next_p++;
	    if(($max_checked + 1) < $next_p){
		$consistent = 0;
		$exact = 0;
	    }
	}
	my $p = $next_p;
	my $match = 0;
	# if its a match only one will match since features can't overlap
	while(($p <= $#$pred_cds) && ($a_stop >= $$pred_cds[$p]->start)){
	    $match = _compare_features($$ann_cds[$a],$$pred_cds[$p]);
	    $nuc_overlap = 1;
	    $max_checked = $p;
	    $p++;
	}
	if($match){
	    $exact_cds = 1;
	}
	else{
	    $consistent = 0;
	    $exact = 0;
	}
    }
    if($max_checked < $#$pred_cds){
	$exact = 0;
	if($ann_high >= 0){
	    $consistent = 0;
	}
    }
    if($exact_cds){
    	$$a_info{Exact_Exon} = 1;
	$$p_info{Exact_Exon} = 1;
    	$$ag_info{Exact_Exon} = 1;
	$$pg_info{Exact_Exon} = 1;
    }
    if($nuc_overlap){
	$$a_info{CDS_Overlap} = 1;
	$$p_info{CDS_Overlap} = 1;
	$$ag_info{CDS_Overlap} = 1;
	$$pg_info{CDS_Overlap} = 1;
    }
    #set consistent and exact values
    if($consistent){
	$$a_info{Consistent} = 1;
	$$p_info{Consistent} = 1;
	$$ag_info{Consistent} = 1;
	$$pg_info{Consistent} = 1;
    }
    if($exact){
	$$a_info{Exact} = 1;
	$$p_info{Exact} = 1;
	$$ag_info{Exact} = 1;
	$$pg_info{Exact} = 1;
	if($utr3_match && $utr5_match){
	    $$a_info{Full_Exact} = 1;
	    $$p_info{Full_Exact} = 1;
	    $$ag_info{Full_Exact} = 1;
	    $$pg_info{Full_Exact} = 1;
	}	
    }
}

sub _compare_features{
    my ($a_feature,$p_feature) = @_;
    my $a_info = $a_feature->tag;
    my $p_info = $p_feature->tag;
    my $a_start = $a_feature->start;
    my $a_stop = $a_feature->stop;
    my $p_start = $p_feature->start;
    my $p_stop = $p_feature->stop;
    my $strand = $a_feature->strand;
    unless($p_feature->strand eq $a_feature->strand){
	return 0;
    }
    if(($a_start == $p_start) &&
       ($a_stop == $p_stop)){
	$$p_info{Correct} = 1; 
	$$p_info{Overlap} = 1;
	$$p_info{Overlap_80p} = 1;
	$$p_info{Splice_5} = 1;
	$$p_info{Splice_3} = 1;
	$$p_info{Correct_Nucs} = [[$p_start,$p_stop]];
	$$a_info{Correct} = 1;
	$$a_info{Overlap} = 1;
	$$a_info{Overlap_80p} = 1;
	$$a_info{Splice_5} = 1;
	$$a_info{Splice_3} = 1;
	$$a_info{Correct_Nucs} = [[$a_start,$a_stop]];
	return 1;
    }
    else{
	my $start = $p_start;
	my $stop = $p_stop;
	if($a_start > $p_start){
	    $start = $a_start;
	}
	if($a_stop < $p_stop){
	    $stop = $a_stop;
	}
	if($start <= $stop){ 
	    if($strand eq "-"){
		if($a_start == $p_start){
		    $$p_info{Splice_3} = 1;
		    $$a_info{Splice_3} = 1;
		}
		if($a_stop == $p_stop){
		    $$p_info{Splice_5} = 1;
		    $$a_info{Splice_5} = 1;
		    
		}
	    }
	    else{
		if($a_start == $p_start){
		    $$p_info{Splice_5} = 1;
		    $$a_info{Splice_5} = 1;
		}
		if($a_stop == $p_stop){
		    $$p_info{Splice_3} = 1;
		    $$a_info{Splice_3} = 1;
		}
	    }
	    my $len = $stop - $start + 1;
	    my $p_len = $p_stop - $p_start + 1;
	    my $a_len = $a_stop - $a_start + 1;
	    $$p_info{Overlap} = 1;
	    $$a_info{Overlap} = 1;
	    if(($len/$p_len >= .8) && ($len/$a_len >= .8)){
		$$p_info{Overlap_80p} = 1;
		$$a_info{Overlap_80p} = 1;
	    }
	    my $a_start = $start;
	    my $a_stop = $stop;
	    my $p_start = $start;
	    my $p_stop = $stop;
	    my $interval = 0;
	    my $new_list = [];
	    while(($interval <= $#{$$p_info{Correct_Nucs}}) &&
		  ($$p_info{Correct_Nucs}[$interval][1] < $p_start)){
		push @$new_list, $$p_info{Correct_Nucs}[$interval];
		$interval++;
	    }
	    while(($interval <= $#{$$p_info{Correct_Nucs}}) &&
		  ($$p_info{Correct_Nucs}[$interval][0] < $p_stop)){
		if($$p_info{Correct_Nucs}[$interval][0] < $p_start){
		    $p_start = $$p_info{Correct_Nucs}[$interval][0];
		}
		if($$p_info{Correct_Nucs}[$interval][1] > $p_stop){
		    $p_stop = $$p_info{Correct_Nucs}[$interval][1];
		}
		$interval++;
	    }
	    push @$new_list, [$p_start,$p_stop];
	    while($interval <= $#{$$p_info{Correct_Nucs}}){
		push @$new_list, $$p_info{Correct_Nucs}[$interval];
		$interval++;
	    }
	    $$p_info{Correct_Nucs} = $new_list;
	    $new_list = [];
	    $interval = 0;
	    while(($interval <= $#{$$a_info{Correct_Nucs}}) &&
		  ($$a_info{Correct_Nucs}[$interval][1] < $a_start)){
		push @$new_list, $$a_info{Correct_Nucs}[$interval];
		$interval++;
	    }
	    while(($interval <= $#{$$a_info{Correct_Nucs}}) &&
		  ($$a_info{Correct_Nucs}[$interval][0] < $a_stop)){
		if($$a_info{Correct_Nucs}[$interval][0] < $a_start){
		    $a_start = $$a_info{Correct_Nucs}[$interval][0];
		}
		if($$a_info{Correct_Nucs}[$interval][1] > $a_stop){
		    $a_stop = $$a_info{Correct_Nucs}[$interval][1];
		}
		$interval++;
	    }
	    push @$new_list, [$a_start,$a_stop];
	    while($interval <= $#{$$a_info{Correct_Nucs}}){
		push @$new_list, $$a_info{Correct_Nucs}[$interval];
		$interval++;
	    }
	    $$a_info{Correct_Nucs} = $new_list;
	}
	return 0;
    }
    return 0;
}

######################################################################
#
#  get_statistics(gtfs,verbose) = @_;
#    Gets general statistics about one or more gtf sets.
#
#    Parameters:
#      gtfs    - reference to a list of GTF objects
#      verbose - boolean turns verbose mode on/off
#
######################################################################
sub get_statistics{
    my ($gtfs,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Statistics Calculations\n";
    }
    my $start_time = time;
    my @data;
    foreach my $gtf_set (@$gtfs){	
	my $gtf_set_genes;
	foreach my $gtf (@$gtf_set){
	    my $gtf_genes = $gtf->genes;
	    push @$gtf_set_genes, $gtf_genes;
	}
	my $data = _get_stats($gtf_set_genes);
	_calculate_stats($data);
	push @data,$data;
    }
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    if($verbose){
	print_time($total_time);
    }
    return @data;
}

# get_statistics internal functions    
sub _get_stats{
    my ($genes) = @_;
    my %data = get_stats_struct();
    foreach my $gene (@$genes){
	_get_gene_list_stats($gene,\%data);
    }
    return \%data;
}

sub _get_gene_list_stats{
    my ($genes,$data) = @_;
    if($#$genes == -1){
	return;
    }
    my $last_stop = $$genes[0]->stop;
    my @pending;
    foreach my $gene (@$genes){
	if($gene->start > $last_stop){
	    _get_gene_stats(\@pending,$data);
	    @pending = ();
	    $last_stop = 0;
	}
	if($gene->stop > $last_stop){
	    $last_stop = $gene->stop;
	}
	push @pending, $gene;
    }
    unless($#pending == -1){
	_get_gene_stats(\@pending,$data);
    }
}

sub _get_gene_stats{
    my ($genes,$data) = @_;
    my @all_txs;
    foreach my $gene (@$genes){
	my @gene_types = get_gene_type_list;
	my $txs = $gene->transcripts;
	my %gene_data = _get_gene_type_hash();
	$gene_data{All} = 1;
	foreach my $type (@gene_types){	    
	    if($gene_data{$type}){
		$$data{Gene}{$type}{Count}++;
		$$data{Gene}{$type}{Total_Transcripts} += $#$txs + 1;
	    }
	}
	push @all_txs, @$txs;
    }
    @all_txs = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_txs;
    _get_tx_stats(\@all_txs,$data);
}

sub _get_tx_stats{
    my ($txs,$data) = @_;
    my @all_cds;
    my @all_utr5;
    my @all_utr3;
    my @all_introns;
    foreach my $tx (@$txs){
	my @tx_types = get_tx_type_list;
	my $cds = $tx->cds;
	my $introns = $tx->introns;
	my $len = $tx->length;
	my $clen = $tx->coding_length;
	my $score = $tx->score;
	my %tx_data = _get_tx_type_hash();
	$tx_data{All} = 1;
	my $starts = $tx->start_codons;
	my $stops = $tx->stop_codons;
	if(($#$starts >= 0) && ($#$stops >= 0)){
	    $tx_data{Complete} = 1;
	    $tx_data{Incomplete} = 0;
	    $$data{Signal}{Start_Codon}{Count}++;
	    $$data{Signal}{Stop_Codon}{Count}++;
	    if($#$cds == 0){
		$tx_data{Single_Exon} = 1;
	    }
	    else{
		$tx_data{Single_Exon} = 0;
	    }		
	}
	else{
	    $tx_data{Complete} = 0;
	    $tx_data{Incomplete} = 1;
	    if($#$starts >=0 ){
		$$data{Signal}{Start_Codon}{Count}++;
	    }
	    elsif($#$stops >= 0){
		$$data{Signal}{Stop_Codon}{Count}++;
	    }
	    else{
		#evan this is an error
	    }
	    $tx_data{Single_Exon} = 0;
	}
	my $total_score = 0;
	foreach my $cds (@$cds){
	    $total_score += $cds->score;
	}
	foreach my $type (@tx_types){
	    if($tx_data{$type}){
		$$data{Transcript}{$type}{Count}++;
		$$data{Transcript}{$type}{Total_Length} += $len;
		push @{$$data{Transcript}{$type}{Length_Array}}, $len;
		$$data{Transcript}{$type}{Total_Coding_Length} += $clen;
		push @{$$data{Transcript}{$type}{Coding_Length_Array}}, $clen;
		$$data{Transcript}{$type}{Total_Exons} += $#$cds + 1;
		push @{$$data{Transcript}{$type}{Exons_Per_Array}}, $#$cds + 1;
		$$data{Transcript}{$type}{Total_Score} += $total_score;
	    }
	}
	push @all_cds, @$cds;
	push @all_utr5, @{$tx->utr5};
	push @all_utr3, @{$tx->utr3};
	push @all_introns, @$introns;
    }
    @all_cds = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_cds;
    @all_utr5 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr5;
    @all_utr3 = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_utr3;
    @all_introns = sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @all_introns;
    _get_exon_stats(\@all_cds,$data);
    _get_exon_stats(\@all_utr5,$data);
    _get_exon_stats(\@all_utr3,$data);
    _get_exon_stats(\@all_introns,$data);
}

sub _get_exon_stats{
    my ($exons,$data) = @_;
    if($#$exons == -1){
	return;
    }
    my @subtypes;
    $subtypes[0]{$$exons[0]->subtype} = 1;
    $subtypes[0]{$$exons[0]->ase} = 1 if ($$exons[0]->ase);
    for(my $i = 0;$i < $#$exons;$i++){
	$subtypes[$i+1]{$$exons[$i+1]->subtype} = 1;
    $subtypes[$i+1]{$$exons[$i+1]->ase} = 1 if ($$exons[$i+1]->ase);
	if($$exons[$i]->equals($$exons[$i+1])){
	    foreach my $subtype (keys %{$subtypes[$i]}){
		$subtypes[$i+1]{$subtype} = 1;
	    }
	    $$exons[$i] = 0;
	}
    }
    my %last_nuc = _get_exon_type_hash();
    for(my $i = 0;$i <= $#$exons;$i++){
	my $exon = $$exons[$i];
	my $subtype_hash = $subtypes[$i];
	unless($exon){
	    next;
	}
	my @exon_types = get_exon_type_list;
	my $len = $exon->length;
	my $score = $exon->score;
	my %exon_data = _get_exon_type_hash();
	foreach my $subtype (keys %$subtype_hash){
	    $exon_data{$subtype} = 1;
	    unless($subtype =~ /Intron/ || $subtype =~ /UTR/){
		$exon_data{All} = 1;
		if(($subtype eq "Initial") ||
		   ($subtype eq "Internal")){
		    $$data{Signal}{Splice_Donor}{Count}++;
		}
		if(($subtype eq "Terminal") ||
		   ($subtype eq "Internal")){
		    $$data{Signal}{Splice_Acceptor}{Count}++;	    
		}
	    }
	}
	foreach my $type (@exon_types){
	    if($exon_data{$type}){
		$$data{Exon}{$type}{Count}++;
		$$data{Exon}{$type}{Total_Length} += $len;
		push @{$$data{Exon}{$type}{Length_Array}}, $len;
		$$data{Exon}{$type}{Total_Score} += $score;
		# count the non-overlapping nucs
		if($exon->stop > $last_nuc{$type}){
		    $$data{Nuc}{$type}{Count} += $len;
		    if($exon->start <= $last_nuc{$type}){
			$$data{Nuc}{$type}{Count} -= $last_nuc{$type} - $exon->start + 1;
		    }
		    $last_nuc{$type} = $exon->stop;
		}
	    }
	}
    }    
}

######################################################################
#
#  get_overlap_statistics($preds,$type,$verbose)
#    Creates clusters of genes, transcrips, or exons which overlap
#    each other by the specified type.  Reports the counts of each cluster
#    type, where a cluster type is a set of clusters containing one or
#    more objects from a particular subset of the preds.  Cluster type
#    counts are returned for all possible subsets of preds.
#
#    Parameters:
#      $preds   - reference to array of sets(arrays) of prediction gtfs
#      $type    - type of overlap.  must come from get_overlap_list()
#      $verbose - boolean turns verbose mode on/off
#
######################################################################
sub get_overlap_statistics{
    my ($preds,$type,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Overlap Calculations\n";
    }
    my $start_time = time;
    my %map = _get_overlap_map();
    my $func = $map{$type};
    my %data =  &$func($preds);
    if($verbose){    
	my $end_time = time;
	my $total_time = $end_time - $start_time + 1;
	print_time($total_time);
    }
    return %data;
}

my @alphabet = qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);

sub _get_overlap_stats{
    my ($preds,$select_func,$compare_func) = @_;
    #create group labels and initialize group data structs
    my @objects;
    my @max;
    for(my $pred = 0;$pred <= $#$preds;$pred++){
	my @pred_objects;
	for(my $gtf_num = 0; $gtf_num <= $#{$$preds[$pred]};$gtf_num++){
	    unless(defined($objects[$gtf_num])){
		$objects[$gtf_num] = [];
	    }
	    my $gtf = $$preds[$pred][$gtf_num];
	    foreach my $object (@{&$select_func($gtf)}){
		$object->set_tag ($alphabet[$pred]);
		push @{$objects[$gtf_num]}, {start => $object->start,
					     stop => $object->stop,
					     list => [$object]};
	    }
	}
    }
    my @clusters;
    my %cluster_count;
    # build the cluster type names
    my @cluster_names;
    for(my $j = 0;$j < 2**($#$preds+1);$j++){
	push @cluster_names, "";
    }
    for(my $i = 0;$i <= $#$preds;$i++){
	my $next = 2**$i;
	my $val = 0;
	for(my $j = 0;$j < 2**($#$preds+1);$j++){
	    if($j == $next){		
		if($val){
		    $val = 0;
		}
		else{
		    $val = 1;
		}
		$next += 2**$i;
	    }
	    if($val){
		$cluster_names[$j] .= $alphabet[$i];
	    }
	}
    }
    foreach my $name (@cluster_names){
	unless($name eq ""){
	    $cluster_count{$name}{total} = 0;
	    for(my $i = 0;$i <= $#$preds;$i++){
		$cluster_count{$name}{$alphabet[$i]} = 0;
	    }
	}
    }
    # build the clusters
    foreach my $set (@objects){
	$set = [sort {$$a{start} <=> $$b{start}} @$set];
	my $next_start = 0;
	for(my $i = 0;$i <= $#$set;$i++){
	    my $object_i = $$set[$i]{list}[0];
	    while(($next_start <= $#$set) &&
		  ($$set[$next_start]{stop} < $$set[$i]{start})){
		_collect_cluster($$set[$next_start],\%cluster_count);
		$next_start++;
	    }
	    my $next = $next_start;
	    while(($next <= $#$set) &&
		  ($$set[$next]{start} <= $$set[$i]{stop})){
		if($next == $i){
		    $next++;
		    next;
		}
		foreach my $object (@{$$set[$next]{list}}){
		    if(&$compare_func($object,$object_i)){
			_combine_clusters($$set[$next],$$set[$i]);
			last;
		    }
		}
		$next++;
	    }
	}
	while($next_start <= $#$set){
	    _collect_cluster($$set[$next_start],\%cluster_count);
	    $next_start++;
	}
    }
    return %cluster_count;
}

sub _collect_cluster{
    my ($cluster,$cluster_count) = @_;
    my @types;
    foreach my $object (@{$$cluster{list}}){
	push @types, $object->tag;
	$object->set_tag();
    }
    @types = sort {$a cmp $b} @types;
    my $last = "";
    my $name = "";
    my %count;
    foreach my $type (@types){
	if(defined($count{$type})){
	    $count{$type}++;
	}
	else{
	    $count{$type} = 1;
	}
	unless($type eq $last){
	    $name .= $type;
	    $last = $type;
	}
    }    
    unless ($name eq ""){
	$$cluster_count{$name}{total}++;
	foreach my $type (keys %count){
	    $$cluster_count{$name}{$type} += $count{$type};
	}
    }
}

sub _combine_clusters{
    my ($c1,$c2) = @_;
    if($$c1{start} < $$c2{start}){
	$$c1{start} = $$c2{start};
    }
    if($$c1{stop} < $$c2{stop}){
	$$c1{stop} = $$c2{stop};
    }
    foreach my $object (@{$$c2{list}}){
	push @{$$c1{list}}, $object;
    }
    $$c2{list} = [];
}

sub get_overlap_labels{
    my ($count) = @_;
    return @alphabet[0..$count];
}

sub get_overlap_mode_list{
    return ("Transcript_Exact_Overlap",
	    "Transcript_Coding_Overlap",
	    "Transcript_Region_Overlap",
	    "Transcript_80p_Region_Overlap",
	    "Transcript_80p_Both_Region_Overlap",
	    "Transcript_Exact_Exon_Overlap",
	    "Transcript_Exact_Intron_Overlap",
	    "Exon_Exact_Overlap",
	    "Exon_One_Base_Overlap",
	    "Exon_80p_Overlap",
	    "Exon_80p_Both_Overlap");
}

sub _get_overlap_map{
    return ("Transcript_Exact_Overlap" => \&get_tx_exact_overlap_statistics,
	    "Transcript_Coding_Overlap" => \&get_tx_coding_overlap_statistics,
	    "Transcript_Region_Overlap" => \&get_tx_1bp_overlap_statistics,
	    "Transcript_80p_Region_Overlap" => \&get_tx_80p_smaller_overlap_statistics,
	    "Transcript_80p_Both_Region_Overlap" => \&get_tx_80p_overlap_statistics,
	    "Transcript_Exact_Exon_Overlap" => \&get_tx_exact_exon_overlap_statistics,
	    "Transcript_Exact_Intron_Overlap" => 
	    \&get_tx_exact_intron_overlap_statistics,
	    "Exon_Exact_Overlap" => \&get_exon_exact_overlap_statistics,
	    "Exon_One_Base_Overlap" => \&get_exon_1bp_overlap_statistics,
	    "Exon_80p_Overlap" => \&get_exon_80p_smaller_overlap_statistics,
	    "Exon_80p_Both_Overlap" => \&get_exon_80p_both_overlap_statistics);
}

sub _get_genes{
    my ($pred) = @_;
    return $pred->genes;
}

sub _get_txs{
    my ($pred) = @_;
    return $pred->transcripts;
}

sub _get_exons{
    my ($pred) = @_;
    my @cds;
    foreach my $tx (@{$pred->transcripts}){
	push @cds, @{$tx->cds};
    }
    return \@cds;
}

sub _1bp_overlap_func{
    my ($a,$b) = @_;
    if(($a->start < $b->stop) &&
       ($a->stop > $b->start) &&
       ($a->strand eq $b->strand)){
	return 1;
    }
    return 0;
}

sub _exact_bounds_overlap_func{
    my ($a,$b) = @_;
    if(($a->start == $b->start) &&
       ($a->stop == $b->stop) &&
       ($a->strand eq $b->strand)){
	return 1;
    }
    return 0;
}

sub _80p_both_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_len = $a->stop - $a->start + 1;
    my $b_len = $b->stop - $b->start + 1;
    my $start = $a->start;
    if($b->start > $start){
	$start = $b->start;
    }
    my $stop = $a->stop;
    if($b->stop < $stop){
	$stop = $b->stop;
    }
    my $overlap = $stop - $start + 1;
    if($overlap > 0){
	if($a_len > $b_len){
	    if(($overlap/$a_len) > .8){
		return 1;
	    }
	}
	else{
	    if(($overlap/$b_len) > .8){
		return 1;
	    }
	}
    }
    return 0;
}

sub _80p_smaller_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_len = $a->stop - $a->start + 1;
    my $b_len = $b->stop - $b->start + 1;
    my $start = $a->start;
    if($b->start > $start){
	$start = $b->start;
    }
    my $stop = $a->stop;
    if($b->stop < $stop){
	$stop = $b->stop;
    }
    my $overlap = $stop - $start + 1;
    if($overlap > 0){
	if($a_len < $b_len){
	    if(($overlap/$a_len) > .8){
		return 1;
	    }
	}
	else{
	    if(($overlap/$b_len) > .8){
		return 1;
	    }
	}
    }
    return 0;
}
sub _tx_coding_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_cds = $a->cds;
    my $b_cds = $b->cds;
    my $i = 0;
    foreach my $cds (@$a_cds){	
	while(($i <= $#$b_cds) && 
	      ($$b_cds[$i]->stop < $cds->start)){
	    $i++;
	}
	if(($i <= $#$b_cds) && 
	   ($$b_cds[$i]->start <= $cds->stop)){	    
	    return 1;
	}
    }
    return 0;
}

sub _tx_exact_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_starts = $a->start_codons;
    my $b_starts = $b->start_codons;
    unless($#$a_starts == $#$b_starts){
	return 0;
    }
    for(my $i = 0;$i <= $#$a_starts;$i++){
	unless(($$a_starts[$i]->start == $$b_starts[$i]->start) &&
	       ($$a_starts[$i]->stop == $$b_starts[$i]->stop)){
	    return 0;
	}
    }
    my $a_stops = $a->stop_codons;
    my $b_stops = $b->stop_codons;
    unless($#$a_stops == $#$b_stops){
	return 0;
    }
    for(my $i = 0;$i <= $#$a_stops;$i++){
	unless(($$a_stops[$i]->start == $$b_stops[$i]->start) &&
	       ($$a_stops[$i]->stop == $$b_stops[$i]->stop)){
	    return 0;
	}
    }
    my $a_exons = $a->cds;
    my $b_exons = $b->cds;
    my $i = 0;
    my $start_i = 0;
    unless($#$a_exons == $#$b_exons){
	return 0;
    }
    for(my $i = 0;$i <= $#$a_exons;$i++){
	unless(($$a_exons[$i]->start == $$b_exons[$i]->start) &&
	       ($$a_exons[$i]->stop == $$b_exons[$i]->stop)){
	    return 0;
	}
    }
    return 1;
}

sub _tx_exact_exon_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_exons = $a->cds;
    my $b_exons = $b->cds;
    my $i = 0;
    my $start_i = 0;
    foreach my $exon (@$a_exons){	
	$i = $start_i;
	while(($i <= $#$b_exons) && 
	      ($$b_exons[$i]->stop < $exon->start)){
	    $i++;
	}
	$start_i = $i;
	while(($i <= $#$b_exons) && 
	      ($$b_exons[$i]->start <= $exon->stop)){	    
	    if(($$b_exons[$i]->start <= $exon->start) &&
	       ($$b_exons[$i]->stop <= $exon->stop)){
		return 1;
	    }
	    $i++;
	}
    }
    return 0;
}

sub _tx_exact_intron_overlap_func{
    my ($a,$b) = @_;
    unless($a->strand eq $b->strand){
	return 0;
    }
    my $a_introns = $a->introns;
    my $b_introns = $b->introns;
    my $i = 0;
    my $start_i = 0;
    foreach my $intron (@$a_introns){	
	$i = $start_i;
	while(($i <= $#$b_introns) && 
	      ($$b_introns[$i]->stop < $intron->start)){
	    $i++;
	}
	$start_i = $i;
	while(($i <= $#$b_introns) && 
	      ($$b_introns[$i]->start <= $intron->stop)){	    
	    if(($$b_introns[$i]->start == $intron->start) &&
	       ($$b_introns[$i]->stop == $intron->stop)){
		return 1;
	    }
	    $i++;
	}
    }
    return 0;
}
 
sub get_tx_exact_exon_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_tx_exact_exon_overlap_func);
}

sub get_tx_exact_intron_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_tx_exact_intron_overlap_func);
}

sub get_tx_coding_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_tx_coding_overlap_func);
}

sub get_tx_1bp_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_1bp_overlap_func);
}

sub get_tx_80p_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_80p_both_overlap_func);
}

sub get_tx_80p_smaller_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_80p_smaller_overlap_func);
}

sub get_tx_exact_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_txs,\&_tx_exact_overlap_func);
}

sub get_exon_1bp_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_exons,\&_1bp_overlap_func);
}

sub get_exon_80p_both_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_exons,\&_80p_both_overlap_func);
}

sub get_exon_80p_smaller_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_exons,\&_80p_smaller_overlap_func);
}

sub get_exon_exact_overlap_statistics{
    my ($preds) = @_;
    return _get_overlap_stats($preds,\&_get_exons,\&_exact_bounds_overlap_func);
}



#evan need to rewrite this
######################################################################
#
#  get_overlap_statisticss($preds,$type,$verbose)
#    Creates clusters of genes, transcrips, or exons which overlap
#    each other by the specified type.  Reports the counts of each cluster
#    which contains only objects from a subset of the predictions
#
#    Parameters:
#      $preds   - reference to array of sets(arrays) of prediction gtfs
#      $type    - type of overlap.  must come from get_overlap_list()
#      $verbose - boolean turns verbose mode on/off
######################################################################
sub filter_predictions{
    my ($ann,$preds,$filter,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Filter Calculations\n";
    }
    my $start_time = time;
    my @filtered_gtfs;
    my @a_genes;
    my @p_genes;
    foreach my $ann_gtf (@$ann){	
	my $genes = $ann_gtf->genes;
	push @a_genes, $genes;
    }
    foreach my $pred (@$preds){
	my $pred_genes;
	foreach my $pred_gtf (@$pred){
	    my $genes = $pred_gtf->genes;
	    push @$pred_genes, $genes;
	}
	push @p_genes, $pred_genes;
	my %stats_struct = get_stats_struct($pred_genes);
    }
    for(my $p = 0;$p <= $#p_genes;$p++){
	my @new_gtfs;
	for(my $i = 0;$i <= $#a_genes;$i++){
	    my $new_genes = [];
	    _filter_gene_lists($a_genes[$i],$p_genes[$p][$i],$new_genes,$filter);
	    my $new_gtf = GTF::new({});
	    foreach my $gene (@$new_genes){
		my $txs = $gene->transcripts;
		my $cds = $gene->cds;
	    }
	    $new_gtf->set_genes($new_genes);
	    push @new_gtfs, $new_gtf;
	}
	push @filtered_gtfs,\@new_gtfs;
    }
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    if($verbose){    
	print_time($total_time);
    }
    return @filtered_gtfs;
}

# functions to compare the data
sub _filter_gene_lists{
    my ($anns,$preds,$new_genes,$filter) = @_;
    my $next_a = 0;
    my $a = 0;
    my $max_init = -1;
    my $last_cleared = -1;
    for(my $p = 0;$p <= $#$preds;$p++){
	_init_gene_tag($$preds[$p]);
	my $p_start = $$preds[$p]->start;
	my $p_stop = $$preds[$p]->stop;
	my $p_strand = $$preds[$p]->strand;	
	while(($next_a <= $#$anns) && ($p_start > $$anns[$next_a]->stop)){
	    if($next_a > $max_init){
		_init_gene_tag($$anns[$next_a]);
		$max_init = $next_a;
	    }
	    $next_a++;
	}
	my $a = $next_a;
	while(($a <= $#$anns) && ($p_stop >= $$anns[$a]->start)){
	    if($a > $max_init){
		_init_gene_tag($$anns[$a]);
		$max_init = $a;
	    }
	    _compare_genes($$anns[$a],$$preds[$p]);
	    $a++;
	}
	my $new_gene = _filter_gene($$preds[$p],$filter);
	if($new_gene){
	    push @$new_genes, $new_gene;
	}
	_clear_gene_tag($$preds[$p]);
	for(my $i = $last_cleared + 1;$i < $next_a;$i++){
	    _clear_gene_tag($$anns[$i]);
	    $last_cleared = $i;
	}
    }   
    for(my $i = $last_cleared + 1;$i <= $max_init;$i++){
	_clear_gene_tag($$anns[$i]);
    }
}

sub _filter_gene{    
    my ($gene,$filter) = @_;
    my $info = $gene->tag;
    my $val = _check_filter($info,"Gene",$filter);
    if($val == -1){
	return 0;
    }
    elsif($val == 0){
	my $new_gene = GTF::Gene::new($gene->id,$gene->seqname,$gene->source,
				      $gene->strand);
	my $txs = $gene->transcripts;
	my $good = 0;
	foreach my $tx (@$txs){
	    my $new_tx = _filter_tx($tx,$filter);
	    if($new_tx){
		$good = 1;
		$new_gene->add_transcript($new_tx);
	    }
	}
	if($good){
	    return $new_gene;
	}
	else{
	    return 0;
	}
    }
    elsif($val == 1){
	my $new_gene = $gene->copy;;
	return $new_gene;
    }
    else{
	print STDERR "Bad val, $val, returned from Eval::_check_filter to ".
	    "Eval::_fitler_gene\n";
	return 0;
    }
}

sub _filter_tx{
    my ($tx,$filter) = @_;
    my $info = $tx->tag;
    my $val = _check_filter($info,"Transcript",$filter);
    if($val == -1){
	return 0;
    }
    elsif($val == 0){
	my $new_tx = GTF::Transcript::new($tx->id);
	my $cds = $tx->cds;
	my $good = 0;
	foreach my $c (@$cds){
	    my $new_cds = _filter_exon($c,$filter);
	    if($new_cds){
		$good = 1;
		$new_tx->add_feature($new_cds);
	    }
	}
	if($good){
	    return $new_tx;
	}
	else{
	    return 0;
	}
    }
    elsif($val == 1){
	my $new_tx = $tx->copy;

	return $new_tx;
    }
    else{
	print STDERR "Bad val, $val, returned from Eval::_check_filter to ".
	    "Eval::_fitler_tx\n";
	return 0;
    }
}

sub _filter_exon{
    my ($exon,$filter) = @_;
    my $info = $exon->tag;
    my $val = _check_filter($info,"Exon",$filter);
    if($val == -1){
	return 0;
    }
    elsif($val == 0){
	return 0;
    }
    elsif($val == 1){
	my $new_exon = $exon->copy();
	return $new_exon;	
    }
    else{
	print STDERR "Bad val, $val, returned from Eval::_check_filter to ".
	    "Eval::_fitler_exon\n";
	return 0;
    }
}

# -1 - fail
#  0 - not yet checked
#  1 - success
#  2 - passed 
sub _check_filter{
    my ($info,$level,$filter) = @_;
    if($$filter[0] eq 'Check'){
	if($level eq $$filter[1]){
	    if($$info{$$filter[2]}){
		return 1;
	    }
	    else{
		return -1;
	    }
	}
	else{
	    if($level eq 'Gene'){
		return 0;
	    }
	    elsif($level eq 'Transcript'){
		if($$filter[1] eq 'Gene'){
		    return 2;
		}
		else{
		    return 0;
		}
	    }
	    else{
		return 2;
	    }
	}
    }
    elsif($$filter[0] eq 'Not'){
	my $val = _check_filter($info,$level,$$filter[1]);
	if($val == 1){
	    return -1;
	}
	elsif($val == -1){
	    return 1;
	}
	return $val;
    }
    elsif($$filter[0] eq 'And'){
	my $val1 = _check_filter($info,$level,$$filter[1]);
	my $val2 = _check_filter($info,$level,$$filter[2]);
	if(($val1 == -1)||($val2 == -1)){
	    return -1;
	}
	elsif(($val1 == 0)||($val2 == 0)){
	    return 0;
	}
	else{
	    return 1;
	}
    }
    elsif($$filter[0] eq 'Or'){
	my $val1 = _check_filter($info,$level,$$filter[1]);
	my $val2 = _check_filter($info,$level,$$filter[2]);
	if(($val1 == 1)||($val2 == 1)){
	    return 1;
	}
	elsif(($val1 == 0)||($val2 == 0)){
	    return 0;
	}
	else{
	    return -1;
	}
    }
    else{
	print STDERR "ERROR in Eval::_check_filter()\n";
	return 0;
    }
}

sub get_filter_types{
    return _get_filter_type_struct();
}

sub _get_filter_type_struct{
    my %filters;
    my @levels = qw(Gene Transcript Exon);
    my @gene_list = get_gene_type_list();
    my @new_list = ();
    foreach my $type (@gene_list){
	if($type ne 'All'){
	    push @new_list, $type;
	}
    }
    @gene_list = @new_list;
    push @gene_list, get_gene_stat_type_list();
    my @tx_list = get_tx_type_list();
    @new_list = ();
    foreach my $type (@tx_list){
	if($type ne 'All'){
	    push @new_list, $type;
	}
    }
    @tx_list = @new_list;
    push @tx_list, get_tx_stat_type_list();
    my @exon_list = get_exon_type_list();
    @new_list = ();
    foreach my $type (@exon_list){
	if($type ne 'All'){
	    push @new_list, $type;
	}
    }
    @exon_list = @new_list;
    push @exon_list, get_exon_stat_type_list();
    return (Levels => \@levels,
	    Gene => \@gene_list,
	    Transcript => \@tx_list,
	    Exon => \@exon_list);
}

######################################################################
#
#  make_graphs($ann,$preds,$x_split_types,$resolution,$verbose)
#    This splits each prediction in $preds according to $resolution 
#    by the $x_split_type and generates eval reports for each split 
#    section of each prediction against the whole annotation.
#
#    Parameters:
#      $ann           - reference to array of annotation gtfs
#      $preds         - reference to array of sets(arrays) of prediction gtfs
#      $graphs        - List of graph types requested.  Graph types are specified 
#                       by a x-split type and a y-level.  Each element of the
#                       list is a hash with two fields {split} and {level} conaining text
#                       values for the x-split and split-level.
#      $reolution     - Hash containing a key for each item in $x_split_type.
#                       Each key points to another hash described below.
#                       hash containing two possible keys.  The first {user}
#                       should be a array of bins for different values of the 
#                       x_split_type.  Each bin is a hash with a {start} and {stop}
#                       field.  The second key is {uniform} and is a hash with 
#                       {min}, {max}, {count}, {size} keys.  Min and max are the 
#                       minimum and maximum values or $x_split_type for any bin,
#                       {count} is the number of bins, and {size} is the size of 
#                       each bin.  All four do not need to be specified.  If any
#                       are skipped they are inferred from the others and the data
#                       if needed.  Only one of the top level keys ({user},{custom})
#                       should be defined.  If both are {user} is used.
#      $verbose       - boolean turns verbose mode on/off
#
######################################################################

sub make_graphs{
    my ($ann,$preds,$graphs,$resolution,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Graph Calculations\n";
    }
    my $start_time = time;
    my %x_splits;
    foreach my $split (get_graph_x_types()){
	$x_splits{$split} = {};
	foreach my $level (get_graph_x_levels()){
	    $x_splits{$split}{$level} = 0;
	}
    }
    foreach my $graph (@$graphs){
	unless(defined($x_splits{$$graph{split}})){
	    print STDERR "Error: Bad x-split value, $$graph{split}, in ".
		"Eval::make_graphs.\n";
	}
	unless(defined($x_splits{$$graph{split}}{$$graph{level}})){
	    print STDERR "Error: Bad split-level value, $$graph{level}, in ".
		"Eval::make_graphs.\n";
	}
	$x_splits{$$graph{split}}{$$graph{level}} = 1;
    }
    my @a_genes;
    my @a_txs;
    my @a_exons;
    my @p_objs;
    my @data;
    foreach my $ann_gtf (@$ann){
	my $genes = $ann_gtf->genes;
	my $txs = $ann_gtf->transcripts;
	my $cds = $ann_gtf->cds;
	push @a_genes, $genes;
	push @a_txs, $txs;
	push @a_exons, $cds;
    }
    my @graph_data;
    # init graph_data struct
    for(my $i = 0;$i <= $#$preds;$i++){
	$graph_data[$i] = {};
	foreach my $split (keys %x_splits){
	    $graph_data[$i]{$split} = {};
	    foreach my $level (keys %{$x_splits{$split}}){
		$graph_data[$i]{$split}{$level} = [];
	    }
	}
    }
    # split preds according to resolution
    my @compare_levels;
    foreach my $split (keys %x_splits){
	my @bins = _get_graph_bins($$resolution{$split},$split);
	foreach my $level (get_graph_x_levels()){
	    if($x_splits{$split}{$level}){
		for(my $i = 0;$i <= $#$preds;$i++){
		    my @split_genes = _split_preds_for_graph($$preds[$i],\@bins,$split,
							     $level);
		    push @p_objs, @split_genes;
		    for(my $bin =0;$bin <= $#split_genes;$bin++){
			my %stats_struct = get_stats_struct();
			push @data, \%stats_struct;
			$graph_data[$i]{$split}{$level}[$bin] = {data => \%stats_struct,
								 min => $bins[$bin]{min},
								 max => 
								     $bins[$bin]{max}};
			push @compare_levels,$level;
		    }
		}
	    }
	}
    }
    # compare gtf sets
    for(my $p = 0;$p <= $#p_objs;$p++){
	for(my $i = 0;$i <= $#a_genes;$i++){
	    if($compare_levels[$p] eq 'Exon'){
		compare_exon_lists($a_exons[$i],$p_objs[$p][$i],$data[$p]);
	    }
	    elsif($compare_levels[$p] eq 'Transcript'){
		compare_tx_lists($a_txs[$i],$p_objs[$p][$i],$data[$p]);
	    }
	    else{
		compare_gene_lists($a_genes[$i],$p_objs[$p][$i],$data[$p]);
	    }
	}
	_calculate_stats($data[$p]);
    }
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    if($verbose){    
	print_time($total_time);
    }
    return @graph_data;    
}
#returning data[pred_number]{x_split}{level}([bin]{~data~})

sub _split_preds_for_graph{
    my ($gtfs,$bins,$split,$level) = @_;    
    my @split;
    foreach my $bin (@$bins){
	my @new;
	foreach my $gtf (@$gtfs){
	    push @new, [];
	}
	push @split,\@new;
    }
    for(my $gtf_num = 0;$gtf_num <= $#$gtfs;$gtf_num++){
	my $gtf = $$gtfs[$gtf_num];
	my @objects;
	if($level eq "Gene"){
	    @objects = @{$gtf->genes};
	}
	elsif($level eq "Transcript"){
	    @objects = @{$gtf->transcripts};
	}
	elsif($level eq "Exon"){
	    @objects = @{$gtf->cds};
	}
	my @vals;
	for(my $i = 0;$i <= $#objects;$i++){
	    my $val = _get_x_val($objects[$i],$split);
	    push @vals, {val => $val,
			 obj => $objects[$i]};
	}
	@vals = sort {$$a{val} <=> $$b{val}} @vals;
	my $bin_num = 0;
	my $val_num = 0;
	while(($val_num <= $#vals) &&
	      $vals[$val_num]{val} < $$bins[$bin_num]{min}){
	    $vals[$val_num]{obj}->set_tag(-1);
	    $val_num++;
	}
	for($val_num = $val_num;$val_num <= $#vals;$val_num++){
	    if($vals[$val_num]{val} <= $$bins[$bin_num]{max}){
		$vals[$val_num]{obj}->set_tag($bin_num);
	    }
	    elsif($bin_num < $#$bins){
		$bin_num++;
		$val_num--;
	    }
	    else{
		last;
	    }
	}
	while($val_num <= $#vals){
	    $vals[$val_num]{obj}->set_tag(-1);
	    $val_num++;
	}
	foreach my $obj (@objects){
	    my $tag = $obj->tag;
	    if($tag != -1){
		push @{$split[$tag][$gtf_num]}, $obj;
	    }
	}
    }
    return @split;
}

sub _get_graph_bins{
    my ($resolution,$split) = @_;
    my @new_bins;
    if(defined($resolution) && defined($$resolution{user})){
	my $bins = $$resolution{user};
	# first check that the bins makes sense.  Fix easily fixable errors
	$bins = [sort {$a->{min} <=> $b->{min}} @$bins];
	my $last_max = -1;
	for(my $i = 0;$i <= $#$bins;$i++){
	    if($$bins[$i]{max} < $$bins[$i]{min}){
		print STDERR "Error: User defined bin has {min} > {max} in for $split ".
		    "in Eval::_get_graph_bins.\n";
		my $temp = $$bins[$i]{min};
		$$bins[$i]{min} = $$bins[$i]{max};
		$$bins[$i]{max} = $temp;
	    }
	    if($i == 0){
		push @new_bins, $$bins[$i];
		$last_max = $$bins[$i]{max};
	    }
	    else{
		if($last_max > $$bins[$i]{min}){
		    print STDERR "Error: Overlapping user defined resolution bins for ".
			"$split in Eval::_get_graph_struct. Bins have been adjusted.\n";
		    $$bins[$i]{min} = $last_max;
		}
		if($last_max < $$bins[$i]{min}){
		    print STDERR "Error: Gaps between user defined resolution bins ".
			"for $split in Eval::_get_graph_struct.  New bin added.\n";
		    push @new_bins, {min => $last_max,
				     max => $$bins[$i]{min}};
		}
		push @new_bins, $$bins[$i];
		$last_max = $$bins[$i]{max};
	    }
	}
    }
    elsif(defined($resolution) && defined($$resolution{uniform})){
	my $data = $$resolution{uniform};
	my $min = 0;
	my $max = 1;
	my $step = .1;
	my $bins = 10;
	my $min_set = 0;
	my $max_set = 0;
	my $step_set = 0;
	# load bins and infer missing values, drop max on overconstrained
	if(defined($$data{min})){
	    $min = $$data{min};
	    $min_set = 1;
	}
	if(defined($$data{max})){
	    $max = $$data{max};
	    if($max < $min){
		my $temp = $min;
		$min = $max;
		$max = $temp;
	    }
	    $step = ($max - $min)/$bins;
	    $max_set = 1;
	}
	if(defined($$data{size})){
	    $step = $$data{size};
	    $bins = ($max - $min)/$step;
	    if($bins > int($bins)){
		$bins++;
	    }
	    $step_set = 0;
	}
	if(defined($$data{count})){
	    $bins = $$data{count};
	    if($max_set){
		if($min_set && !$step_set){
		    $step = ($max - $min)/$bins;
		}
		elsif(!$min_set && $step_set){
		    $min = $max - $bins*$step;		    
		}
	    }
	}
	$max = $min + $step;
	for(my $bin = 0;$bin <= $bins;$bin++){
	    $new_bins[$bin] = {min => $min,
			       max => $max};
	    $min += $step;
	    $max += $step;
	}
    }
    else{
	print STDERR "Error: Missing resolution for $split in Eval::".
	    "_get_graph_struct.\n";
	my $bins = 10;
	my $min = 0;
	my $step = 1000;
	my $max = $min + $step;
	for(my $bin = 0;$bin <= $bins;$bin++){
	    $new_bins[$bin] = {min => $min,
			       max => $max};
	    $min += $step;
	    $max += $step;
	}    
    }
    return @new_bins;
}

sub _get_x_val{
    my ($obj,$type) = @_;
    my %map = _get_graph_x_val_map();
    my $func = $map{$type};
    return &$func($obj);
}

sub get_graph_x_types{
    my %map = _get_graph_x_val_map();
    my @list;
    foreach my $type (keys %map){
	push @list, $type;
    }
    return @list;
}

sub _get_graph_x_val_map{
    return ("GC%" => \&_get_gc_percent,
	    "Match%" => \&_get_match_percent,
	    "Mismatch%" => \&_get_mismatch_percent,
	    "Unaligned%" => \&_get_unaligned_percent,
	    "Length" => \&_get_length);
}

sub _get_gc_percent{
    my ($obj) = @_;
    return $obj->gc_percentage;
}

sub _get_match_percent{
    my ($obj) = @_;
    return $obj->match_percentage;
}

sub _get_mismatch_percent{
    my ($obj) = @_;
    return $obj->mismatch_percentage;
}

sub _get_unaligned_percent{
    my ($obj) = @_;
    return $obj->unaligned_percentage;
}

sub _get_length{
    my ($obj) = @_;
    return $obj->length;
}

sub get_graph_y_types{
    return get_list_struct();
}

sub get_graph_x_levels{   
    return qw(Gene Transcript Exon);
}


######################################################################
#
#  get_distribution($gtf,$distribution)
#    This function calculates each distribution in the $distributions
#    list for each gtf set in the $gtfs list; 
#
#    Parameters:
#      $gtfs          - A reference to an array of sets of gtfs.  Each 
#                       gtf set is an array of gtfs.
#      $distributions - A refence to an array of distribution types.
#                       Each distribution type is a string form the list
#                       Returned by get_distribution_type_list;
#      $verbose       - boolean turns verbose mode on/off
#
######################################################################

sub get_distribution{
    my ($gtfs,$distributions,$verbose) = @_;
    unless(defined($verbose)){
	$verbose = 0;
    }
    if($verbose){
	print STDERR "Starting Filter Calculations\n";
    }
    my $start_time = time;
    my %gene_funcs;
    my %tx_funcs;
    my %exon_funcs;
    my $dist_funcs = _get_distribution_functions();
    my @data;
    for(my $i = 0;$i <= $#$gtfs;$i++){
	$data[$i] = {};
	foreach my $dist (keys %$distributions){
	    if($$distributions{$dist} != 0){
		$data[$i]{$dist} = {};
	    }
	}
    }
    foreach my $dist (keys %$distributions){
	if($$distributions{$dist} != 0){
	    if($dist =~ /Gene/){
		$gene_funcs{$dist} = $$dist_funcs{$dist};
	    }
	    elsif($dist =~ /Transcript/){
		$tx_funcs{$dist} = $$dist_funcs{$dist};
	    }
	    elsif($dist =~ /Exon/){
		$exon_funcs{$dist} = $$dist_funcs{$dist};
	    }
	}
    }
    my @size = keys %gene_funcs;
    if($#size >= 0){
	_get_distribution($gtfs,\@data,\&_get_genes,\%gene_funcs);
    }
    @size = keys %tx_funcs;
    if($#size >= 0){
	_get_distribution($gtfs,\@data,\&_get_txs,\%tx_funcs);	
    }
    @size = keys %exon_funcs;
    if($#size >= 0){
	_get_distribution($gtfs,\@data,\&_get_exons,\%exon_funcs);	
    }
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    if($verbose){    
	print_time($total_time);
    }
    return @data;
}

sub _get_distribution{
    my ($gtfs,$data,$type_func,$dist_funcs) = @_;
    for(my $i = 0;$i <= $#$gtfs;$i++){
	my $list = $$gtfs[$i];
	my @objects;
	foreach my $gtf (@$list){
	    my $list_objects = &$type_func($gtf);
	    push @objects, @$list_objects;
	}
	foreach my $object (@objects){
	    foreach my $dist_name (keys %$dist_funcs){
		my $func = $$dist_funcs{$dist_name};
		my $val = &$func($object);
		if(defined($$data[$i]{$dist_name}{$val})){
		    $$data[$i]{$dist_name}{$val}++;
		}
		else{
		    $$data[$i]{$dist_name}{$val} = 1;
		}
	    }
	}
    }
}

sub get_distribution_type_hash{
    my $dists = _get_distribution_functions();
    foreach my $dist (keys %$dists){
	$$dists{$dist} = 0;
    }
    return %$dists;
}

sub get_distribution_type_list{
    return qw(Transcripts_Per_Gene
	      Transcript_Length
	      Transcript_Coding_Length
	      Exons_Per_Transcript
	      Exon_Length
	      Exon_Score);
}

sub _get_distribution_functions{
    return ({Transcripts_Per_Gene => \&_get_txs_per,
	     Transcript_Length => \&_get_tx_length,
	     Transcript_Coding_Length => \&_get_coding_length,
	     Exons_Per_Transcript => \&_get_exons_per,
	     Exon_Length => \&_get_exon_length,
	     Exon_Score => \&_get_exon_score});
}

sub _get_txs_per{
    my ($gene) = @_;
    my $txs = $gene->transcripts;
    return($#$txs+1);
}

sub _get_exons_per{
    my ($tx) = @_;
    my $exons = $tx->cds;
    return($#$exons+1);
}

sub _get_tx_length{
    my ($tx) = @_;
    return($tx->coding_stop - $tx->coding_start +1);
}

sub _get_coding_length{
    my ($tx) = @_;
    return($tx->coding_length);
}

sub _get_exon_length{
    my ($exon) = @_;
    return($exon->stop - $exon->start + 1);
}

sub _get_exon_score{
    my ($exon) = @_;
    return($exon->score);
}

######################################################################
#
#  General Functions
#
######################################################################

sub print_time{
    my ($total_time) = @_;
    my ($days,$hours,$minutes,$seconds) = (0,0,0,0);
    if($total_time > 86400){
	$days = int($total_time/86400);
	$total_time = $total_time % 86400;
    }
    if($total_time > 3600){
	$hours = int($total_time/3600);
	$total_time = $total_time % 3600;
    }
    if($total_time > 60){
	$minutes = int($total_time/60);
	$total_time = $total_time % 60;
    }
    $seconds = $total_time;    
    print STDERR "Calculations Complete\n";
    print STDERR "\tTime:";
    if($days > 0){
	print STDERR " $days Days,";
    }
    if($hours > 0){
	print STDERR " $hours Hours,";
    }
    if($minutes > 0){
	print STDERR " $minutes Minutes,";
    }
    print STDERR " $seconds Seconds\n";
}

sub max{
    my ($a,$b) = @_;
    if($a > $b){
	return $a;
    }
    return $b;
}

# sort functions (faster?)

sub sort{
    my ($decide_func,@data) = @_;
    print "running my sort\n";
    my @temp = @data;
    qSort2($decide_func,\@temp);
    return @temp;
}

sub qSortPart{
    my ($decide_func,$lo,$hi,$data) = @_;
    if($lo >= $hi){
	return;
    }
    if($lo + 1 == $hi){
	if(&$decide_func($$data[$lo],$$data[$hi]) == -1){
	    my $temp = $$data[$lo];
	    $$data[$lo] = $$data[$hi];
	    $$data[$hi] = $temp;
	}
	return;
    }
    my $tLo = $lo;
    my $tHi = $hi;
    my $tMid = int(($lo+$hi)/2);
    my $pivot = $$data[$tMid];
    $$data[$tMid] = $$data[$hi];
    $$data[$hi] = $pivot;
    while($tLo < $tHi){
	while((&$decide_func($$data[$tLo],$pivot) != -1) && ($tLo < $tHi)){
	    $tLo++;
	}
	while((&$decide_func($$data[$tHi],$pivot) != 1) && ($tLo < $tHi)){
	    $tHi--;
	}
	if($tLo < $tHi){
	    my $temp = $$data[$tLo];
	    $$data[$tLo] = $$data[$tHi];
	    $$data[$tHi] = $temp;
	}
    }
    
    $$data[$hi] = $$data[$tHi];
    $$data[$tHi] = $pivot;
    qSortPart($decide_func,$lo,$tLo-1,$data);
    qSortPart($decide_func,$tHi+1,$hi,$data);
}


sub qSort2{
    my ($decide_func,$data) = @_;
    my $q = 1;
    my $queue = [[0,$#$data]];
    while($q > 0){
	my $lo = $$queue[$q-1][0];
	my $hi = $$queue[$q-1][1];
	$q--;
	if($lo >= $hi){
	    next;
	}
	if($lo + 1 == $hi){
	    if(&$decide_func($$data[$lo],$$data[$hi]) == -1){
		my $temp = $$data[$lo];
		$$data[$lo] = $$data[$hi];
		$$data[$hi] = $temp;
	    }
	    next;
	}
	my $tLo = $lo;
	my $tHi = $hi;
	my $tMid = int(($lo+$hi)/2);


	my $pivot = $$data[$tMid];
	$$data[$tMid] = $$data[$hi];
	$$data[$hi] = $pivot;
	while($tLo < $tHi){
	    while((&$decide_func($$data[$tLo],$pivot) != -1) && ($tLo < $tHi)){
		$tLo++;
	    }
	    while((&$decide_func($$data[$tHi],$pivot) != 1) && ($tLo < $tHi)){
		$tHi--;
	    }
	    if($tLo < $tHi){
		my $temp = $$data[$tLo];
		$$data[$tLo] = $$data[$tHi];
		$$data[$tHi] = $temp;
	    }
	}
	
	$$data[$hi] = $$data[$tHi];
	$$data[$tHi] = $pivot;
	
	$$queue[$q] = [$lo, $tLo -1];
	$q++;	
	$$queue[$q] = [$tHi+1,$hi];
	$q++;
    }
}

1;
__END__


sub count_loci{
    my ($ann) = @_;
    my %overlap;
    my $loci_count = 0;
    my @strands = ('+','-');
    foreach my $gtf (@$ann){
	my $genes = $gtf->genes;
	foreach my $strand (@strands){
	    for(my $i = 0;$i <= $#$genes;$i++){
		unless($$genes[$i]->strand eq $strand){
		    next;
		}
		$loci_count++;
		my $j = $i + 1;
		while(($j <= $#$genes) &&
		      ($$genes[$i]->stop >= $$genes[$j]->start)){
		    $j++;
		}
		$i = $j-1;
	    }
	}
    }
    return $loci_count;
}
