#!/usr/local/bin/perl

#Copyright (c) 2003 by Mihaela Pertea.

use strict;
use FileHandle;

my (%vector,@ag,@gt,@scoreag,@scoregt);
my (@gts,@ags,@sgt,@sag);    
my (@neg,@pos,@score);

my $no_leaves=10;
my $no_of_trees=10;

return 1;

# this procedure computes the *.false file by evaluating scores
# based on their rank in the false.* file
sub calcfalse {
    my ($input,$output,$ind)=@_;

    open(I,$input);
    my $i=0;
    my @score;
    my @neg;
    my @poz;
    my @negperc;
    my @pozperc;
    my ($totalpos,$totalneg);

    <I>;

    while(<I>) {
	chomp;
	my @a=split;
	$score[$i]=$a[0];
	$neg[$i]=$a[1];
	chop($a[3]);
	$negperc[$i]=$a[3];
	$poz[$i]=$a[5];
	chop($a[7]);
	if(!$totalpos && $a[7]!=0.0) { $totalpos=100*$a[5]/$a[7];}
	$pozperc[$i++]=$a[7];
    }
    
    $totalneg=$neg[$#score];
    close(I);

    open(O,">$output");

    print O $#score+1,"\n";
    $poz[$#score+1]=0;
    for($i=0;$i<=$#score;$i++) {
	if($ind && $totalpos) { # var for acceptor and donors
	    print O $score[$i]," ",($i+1)*(1-$poz[$i]/$totalpos)/($#score+2),"\n";
	}
	else { # var for start and stop sites
	    print O $score[$i]," ",($i+1)/($#score+2),"\n"; 
	}
    }

    close(O);
}



# this procedure takes a false.* file and transfomrs it into a DT file
sub false2dt {
    my ($input,$output)=@_;

    my $totalneg=0;
    my $totalpos=0;
    my $prevno=0;
    my $i=0;
    my ($lastnegno,$lastnegperc);

    open(I,$input);

    <I>;
    while(<I>) {
	chomp;
	my @a=split;
	if(!$totalpos) { 
	    chop($a[7]);
	    $totalpos=int($a[5]*100/$a[7]);
	}
	if($a[5]!=$prevno){
	    $score[$i]=$a[0];
	    $neg[$i]=$a[1];
	    $pos[$i++]=$a[5];
	    $prevno=$a[5];
	    chop($a[3]);
	    $lastnegperc=$a[3];
	    $lastnegno=$a[1];
	}
    }

    close(I);

    $totalneg=int($lastnegno*100/$lastnegperc);

    open(O,">$output");
    print O "Training set: file.train, Dimensions: 1, Categories: 2\n\n\n";
    my $totno=$i;

#err
#for($i=0;$i<$totno;$i++) {
#    print "$i ",$score[$i]," ",$neg[$i]," ",$pos[$i],"\n";
#}


    $i=int($totno/2);
    my $leftnegtotal=$neg[$i];
    my $leftpostotal=$totalpos-$pos[$i];
    my $rightnegtotal=$totalneg-$neg[$i];
    my $rightpostotal=$pos[$i];
    print O "Root Hyperplane: Left = [$leftnegtotal,$leftpostotal], Right = [$rightnegtotal,$rightpostotal]\n";
    print O "1.000000 x[1] + ",-$score[$i]," = 0\n\n";

    # left case                                                                                                                 
    if($i-1>=0 && $leftnegtotal && $leftpostotal) {
	my $node="l";
	formtree(*O,$node,0,$i-1,0,$totalpos,$neg[$i],$pos[$i]);
    }
    # right case                                                                                                                
    if($i+1<=$totno-1 && $rightnegtotal && $rightpostotal) {
	my $node="r";
	formtree(*O,$node,$i+1,$totno-1,$neg[$i],$pos[$i],$totalneg,0);
    }
    close(O);
}

# this procedure is called by false2dt to form a subtree from
# a given node
sub formtree {
    my ($out,$node,$left,$right,$negleft,$posleft,$negright,$posright)=@_;

    my $i=$left+int(($right-$left+1)/2);

    #err
    #print "midpoint =$i\n";

    my $leftnegtotal=$neg[$i]-$negleft;
    my $leftpostotal=$posleft-$pos[$i];
    my $rightnegtotal=$negright-$neg[$i];
    my $rightpostotal=$pos[$i]-$posright;

    print $out "$node Hyperplane: Left = [$leftnegtotal,$leftpostotal], Right = [$rightnegtotal,$rightpostotal]\n";
    print $out "1.000000 x[1] + ",-$score[$i]," = 0\n\n";

    if($left!=$right) {
	# left case
	if($i-1>=$left && $leftnegtotal && $leftpostotal) {
	    my $tempnode=$node."l";

	    #err
#	    print "call to $tempnode ",$left," ",$i-1," ",$negleft," ",$negright," ",$neg[$i]," ",$pos[$i],"\n";

	    formtree($out,$tempnode,$left,$i-1,$negleft,$posleft,$neg[$i],$pos[$i]);
	}
	# right case
	if($i+1<=$right && $rightnegtotal && $rightpostotal) {
	    my $tempnode=$node."r";

	    #err
#	    print "call to $tempnode ",$i+1," ",$right," ",$negleft," ",$negright," ",$neg[$i]," ",$pos[$i],"\n";

	    formtree($out,$tempnode,$i+1,$right,$neg[$i],$pos[$i],$negright,$posright);
	}
    }
}



# this procedure takes exons.dat, genelist and seqs as arguments
# and computes the in-frame hexamer frequencies for exons and the 
# hexamer frequencies for all frames of exons, introns, and intergenic
# DNA; the results are written in 2 files: exon.hexfreq and 
# all.hexfreq
sub hexfreq {

    my ($f,$g,$h)=@_;

    my ($allfreq,$exonfreq,%takegene,%exons);

    # read genes to be considered
    open(F,$g) || die "ERROR 15: Could not open file $g for reading!\n";

    while(<F>) {
	chomp;
	$takegene{$_}++;
    }
    close(F);

    # read exons
    open(F,$f) || die "ERROR 16: Could not open file $f for reading!\n";

    while(<F>) {
	chomp;
	if($_) {
	    my @a=split;
	    if($takegene{$a[0]}) {
		$exons{$a[0]}.=$a[1]." ".$a[2]." ";
	    }
	}
    }

    close(F);

    # initialize hexamer frequencies

    my (@let,$hexamer,%afreq,%efreq);
    $let[0]='A';
    $let[1]='C';
    $let[2]='G';
    $let[3]='T';

    for(my $i1=0;$i1<4;$i1++) {
	for(my $i2=0;$i2<4;$i2++) {
	    for(my $i3=0;$i3<4;$i3++) {
		for(my $i4=0;$i4<4;$i4++) {
		    for(my $i5=0;$i5<4;$i5++) {
			for(my $i6=0;$i6<4;$i6++) {
			    $hexamer=$let[$i1].$let[$i2].$let[$i3].$let[$i4].$let[$i5].$let[$i6];
			    $afreq{$hexamer}=0;
			    $efreq{$hexamer}=0;
			}
		    }
		}
	    }
	}
    }
    

    # compute hexamer frequencies

    open(F,$h) || die "ERROR 17: Could not open file $h for reading!\n";

    while(<F>){
	chomp;
	my @a=split;

	if($takegene{$a[0]}) {
	    my $len=length($a[1]);
	    for(my $i=0;$i+5<$len;$i++) {
		$hexamer=uc(substr($a[1],$i,6));

		$afreq{$hexamer}++;
	    }
	
	    my @valex=split(/ /,$exons{$a[0]});
	    my $frame=0;

	    my ($exon);

	    for(my $i=0;$i<=$#valex;$i+=2) {
		$exon=substr($a[1],$valex[$i]-1,$valex[$i+1]-$valex[$i]+1);
		$len=length($exon);
		for(my $j=$frame;$j+5<$len;$j+=3) { # here I compute the in-frame frequency
		    $hexamer=uc(substr($exon,$j,6));
		    $efreq{$hexamer}++;
		}
		$len-=$frame;
		$frame=(3-($len%3))%3;
	    }
	}
    }

    close(F);

    my $exonfile="exon.hexfreq";
    my $allfile="all.hexfreq";

    for(my $i1=0;$i1<4;$i1++) {
	for(my $i2=0;$i2<4;$i2++) {
	    for(my $i3=0;$i3<4;$i3++) {
		for(my $i4=0;$i4<4;$i4++) {
		    for(my $i5=0;$i5<4;$i5++) {
			for(my $i6=0;$i6<4;$i6++) {
			    
			    $hexamer=$let[$i1].$let[$i2].$let[$i3].$let[$i4].$let[$i5].$let[$i6];
			    $allfreq+=$afreq{$hexamer};
			    $exonfreq+=$efreq{$hexamer};
			}
		    }
		}
	    }
	}
    }

    open(E,">$exonfile") or die "ERROR 18: Could not open file $exonfile for writing!\n";
    open(A,">$allfile") or die "ERROR 19: Could not open file $allfile for writing!\n";

    

    for(my $i1=0;$i1<4;$i1++) {
	for(my $i2=0;$i2<4;$i2++) {
	    for(my $i3=0;$i3<4;$i3++) {
		for(my $i4=0;$i4<4;$i4++) {
		    for(my $i5=0;$i5<4;$i5++) {
			for(my $i6=0;$i6<4;$i6++) {
			    $hexamer=$let[$i1].$let[$i2].$let[$i3].$let[$i4].$let[$i5].$let[$i6];
			    my $freq=$efreq{$hexamer}/$exonfreq;
			    printf E "%s %.10f\n",$hexamer,$freq;
			    $freq=$afreq{$hexamer}/$allfreq;
			    printf A "%s %.10f\n",$hexamer,$freq;
			}
		    }
		}
	    }
	}
    }
    close(E);
    close(A);
}


# this procedure takes exon.hexfreq, all.hexfreq (the files generated by 
# by hexfreq), genelist, exons.dat and sites.scores (the file 
# generated by getsites.pl) as arguments and extract the features for
# the true and false exons; the results are written in 4 files, 
# according with the exon types: firstexons.features, snglexons.features
# internalexons.features and lastexons.features
# as well as the intron features written in introns.features
sub feature_extract {

    my ($c,$d,$e,$f,$g,$h,$max_ex,$max_in)=@_;

    # read genes to be considered

    my %takegene;

    open(F,$c) || die "ERROR 1: Could not open file $c for reading: $!\n";    # genelist
    while(<F>) {
	chomp;
	$takegene{$_}++;
    }
    close(F);

    # read frequencies

    my (%ehexamer,%ahexamer);

    open(F,$d) || die "ERROR 2: Could not open file $d for reading: $!\n";    # exon.hexfreq
    while(<F>) {
	chomp;
	if($_) {
	    my @a=split;
	    $ehexamer{$a[0]}=$a[1];
	}
    }
    close(F);

    open(F,$e) || die "ERROR 3: Could not open file $e for reading: $!\n";    # all.hexfreq

    while(<F>) {
	chomp;
	if($_) {
	    my @a=split;
	    $ahexamer{$a[0]}=$a[1];
	}
    }
    close(F);

    my (@let);

    $let[0]='A';
    $let[1]='C';
    $let[2]='G';
    $let[3]='T';

    for(my $i1=0;$i1<4;$i1++) {
	for(my $i2=0;$i2<4;$i2++) {
	    for(my $i3=0;$i3<4;$i3++) {
		for(my $i4=0;$i4<4;$i4++) {
		    for(my $i5=0;$i5<4;$i5++) {
			for(my $i6=0;$i6<4;$i6++) {
			    my $hexamer=$let[$i1].$let[$i2].$let[$i3].$let[$i4].$let[$i5].$let[$i6];
			    if($ehexamer{$hexamer}==0) {
				$vector{$hexamer}=-15;}
			    else {
				$vector{$hexamer}=log($ehexamer{$hexamer}/$ahexamer{$hexamer});
			    }
			}
		    }
		}
	    }
	}
    }

    # read exons

    my %exons;

    open(F,$f) || die "ERROR 4: Could not open file $f for reading: $!\n";    # exons.dat
    while(<F>) {
	chomp;
	if($_) {
	    my @a=split;
	    if($takegene{$a[0]}) {
		$exons{$a[0]}.=" ".$a[1]." ".$a[2]." ";
	    }
	}
    }
    close(F);


    # read sequences

    my %seqs;

    open(F,$g) || die "ERROR 5: Could not open file $g for reading: $!\n";    # seqs
    while(<F>) {
	chomp;
	my ($name,$seq) = /^(\S+)\s+(\S+)$/;
	$seqs{$name}=$seq;
    }
    close(F);
    
#     open(F,$w) ||  die "ERROR 169: Could not open file $satg for reading: $!\n"; # sites.atg
#     <F>;
#     my $true=1;
#     while(<F>) {
# 	chomp;
# 	if($_) {
# 	    if($_ eq "Scores for false starts") { $true=0; }
# 	    else {
# 		my @a=split;
# 		if($true) {
# 		    $truestart{$a[0]}=$a[1]." ".$a[2];
# 		}
# 		else { 
# 		    $falsestart{$a[0]}.=$a[1]." ".$a[2]." ";
# 		}
# 	    }
# 	}
#     }
#     close(F);
# 
# 
#     open(F,$z) ||  die "ERROR 170: Could not open file $sstop for reading: $!\n"; # sites.stop
#     <F>;
#     $true=1;
#     while(<F>) {
# 	chomp;
# 	if($_) {
# 	    if($_ eq "Scores for false stops") { $true=0; }
# 	    else {
# 		my @a=split;
# 		if($true) {
# 		    $truestop{$a[0]}=$a[1]." ".$a[2];
# 		}
# 		else { 
# 		    $falsestop{$a[0]}.=$a[1]." ".$a[2]." ";
# 		}
# 	    }
# 	}
#     }
#     close(F);

    # print the features

#     my $fe="firstexons.features";
    my $ie="internalexons.features";
#     my $le="lastexons.features";
#     my $intr="introns.features";
    my $feall="firstexons.features.allinfo";
    my $leall="lastexons.features.allinfo";

#     open(FEx,">$fe");
    open(IEx,">$ie");
#     open(LEx,">$le");
#     open(Intr,">$intr");
    open(FExAll,">$feall");
    open(LExAll,">$leall");

    open(F,$h) || die "ERROR 6: Could not open file $h for reading: $!\n";     # sites.scores

    my ($name,$lastname,$iag,$igt,$itag,$itgt,$ifag,$ifgt);

    $lastname="";

    while(<F>) {
	chomp;
	if(($name) = /^Gene (\S+)$/) {
	    if($takegene{$lastname}) {
		print "Process $lastname...\n";
		processgene($lastname,$exons{$lastname},$iag,$igt,$seqs{$lastname},$max_ex,$max_in);		
	    }
	
	    $iag=-1;$itag=0;$ifag=0;
	    $igt=-1;$itgt=0;$ifgt=0;
	    if($takegene{$name}) { $lastname=$name; }
	    else { $lastname=""; }
	}
	elsif($lastname) {

	    my @a=split;
	    if($a[0] eq "ag:") {
		$name=$a[1]+2;
		$name=" ".$name." ";
		if(index($exons{$lastname},$name)>-1) {
		    $ag[++$iag]=$a[1];
		    $scoreag[$iag]=$a[2];
		    $itag++;
		}
		else {
		    my $yes=1;
		    if($ifag>$itag) {
			my $roll = int(rand 20);
			$yes=0 if $roll;
		    }
		    if($yes) {
			$ag[++$iag]=$a[1];
			$scoreag[$iag]=$a[2];
			$ifag++;
		    }
		}
	    }
	    if($a[0] eq "gt:") {
		$name=$a[1]-1;
		$name=" ".$name." ";
		if(index($exons{$lastname},$name)>-1) {
		    $gt[++$igt]=$a[1];
		    $scoregt[$igt]=$a[2];
		    $itgt++;
		}
		else {
		    my $yes=1;
		    if($ifgt>$itgt) {
			my $roll = int(rand 20);
			$yes=0 if $roll;
		    }
		    if($yes) {
			$gt[++$igt]=$a[1];
			$scoregt[$igt]=$a[2];
			$ifgt++;
		    }
		}
	    }
	}

    }


    if($takegene{$lastname}) {
	print "Process $lastname...\n";
	processgene($lastname,$exons{$lastname},$iag,$igt,$seqs{$lastname},$max_ex,$max_in);		
    }
    
    close(F);
#    FEx->autoflush(1);
#    close(FEx);
    close(FExAll);
    IEx->autoflush(1);
    close(IEx);
#    LEx->autoflush(1);
#    close(LEx);
    close(LExAll);
#    Intr->autoflush(1);
#    close(Intr);

    # process single genes

#    my $singlex="snglexons.features";
    my $singlexall="snglexons.features.allinfo";
#    open(Ex,">$singlex");
    open(ExAll,">$singlexall");

    foreach $name (keys %takegene) {
	print "Process $name...\n";
	
	my @ex=split(/\s+/,$exons{$name});

	if($#ex==2) {
# 	    my $length=$ex[2]-$ex[1]+1;
# 	    my $exon=substr($seqs{$name},$ex[1]-1,$length-3);
# 	    my $hexv=hexval($exon);
# 	    print Ex "$length $hexv 1\n";
	}
	else {
	    my $where=$ex[1]-1;
	    my $cod=substr($seqs{$name},$where,3);
	    my $exon="";
	    while($cod && ($cod ne "tga")&&($cod ne "taa")&&($cod ne "tag")) {
		$exon.=$cod;
		$where+=3;
		$cod=substr($seqs{$name},$where,3);
	    }
	    my $length=length($exon);
	    if($length>=60) {
		my $hexv=hexval($exon);
#		print Ex "$length $hexv 2\n";
		print ExAll "$name $length $hexv ",$ex[1]," 2\n";
	    }
	}
    }

#     Ex->autoflush(1);
#     close(Ex);

}


# oneorf is used by feature_extract 
sub oneorf {
    my $prelex=$_[0];
    my $lengthex=length($prelex);
    
    my $ret=3;

    for(my $frame=0;$frame<3;$frame++) {
	for(my $pos=$frame;$pos+2<$lengthex;$pos+=3)  {
	    my $cod=substr($prelex,$pos,3);
	    if(($cod eq "taa")||($cod eq "tag")||($cod eq "tga")) {
		$ret--;
		last;
	    }
	}
    }

    return($ret);
}

# hexval is used by feature_extract
sub hexval {
    my $prelex=$_[0];
    my $lengthex=length($prelex);

    my @hex;

    for(my $frame=0;$frame<3;$frame++) {
	$hex[$frame]=0;
	for(my $pos=$frame;$pos+5<$lengthex;$pos+=3)  {
	    my $hexamer=uc(substr($prelex,$pos,6));
	    $hex[$frame]+=$vector{$hexamer};
	}
    }

    my $max=$hex[0];
    if($hex[1]>$max) { $max=$hex[1]; }
    if($hex[2]>$max) { $max=$hex[2]; }
	    
    return($max);
}

sub processgene {
   
    my ($name,$exons,$iag,$igt,$seq,$max_ex,$max_in)=@_;
 
    # process gene
    
    my @ex=split(/\s+/,$exons);
    
    if($#ex>2) {
	my $firstgt=$ex[2]+1;
	my $lastag=$ex[$#ex-1]-2;
	my $firstexon=1;

	my ($orf0,$orf1,$orf2,$process,$lenprocess,$last,$length,$exon,$intron,$hexv);
		
	for(my $i=0;$i<=$iag;$i++) {
	    
	    # check here for last exons
	    if($lastag != $ag[$i]) {
		$orf0=1;
		$orf1=1;
		$orf2=1;
		$process=substr($seq,$ag[$i]+1);
		$lenprocess=length($process);
		$last=$ag[$i]+1;
		
		while($orf0 || $orf1 || $orf2) {
		    
		    $length=$last+2-$ag[$i];
		    if($length > $max_ex || $length+3>$lenprocess) { last; }
		    
		    if($orf0) {
			my $cod0=substr($process,0,3);
			if(($cod0 eq "taa")||($cod0 eq "tag")||($cod0 eq "tga")) { 
			    $orf0=0;
			    if($length-3<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$ag[$i]+1,$length-3);
				$hexv=hexval($exon);
			    }
#			    print LEx $scoreag[$i]," ",$length-3," ",$hexv," 2\n";
			    print LExAll $name," ",$scoreag[$i]," ",$length-3," ",$hexv," ",$ag[$i]+2+$length-3," 2\n";
			}
		    }
		
		    if($orf1) {
			my $cod1=substr($process,1,3);
			$length=$last+2-$ag[$i]+1;
			if(($cod1 eq "taa")||($cod1 eq "tag")||($cod1 eq "tga")) { 
			    $orf1=0;
			    if($length-3<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$ag[$i]+1,$length-3);
				$hexv=hexval($exon);
			    }
# 			    print LEx $scoreag[$i]," ",$length-3," ",$hexv," 2\n";
			    print LExAll $name," ",$scoreag[$i]," ",$length-3," ",$hexv," ",$ag[$i]+2+$length-3," 2\n";
			}
		    }
		
		    if($orf2) {
			my $cod2=substr($process,2,3);
			$length=$last+2-$ag[$i]+2;
			if(($cod2 eq "taa")||($cod2 eq "tag")||($cod2 eq "tga")) { 
			    $orf2=0;
			    if($length-3<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$ag[$i]+1,$length-3);
				$hexv=hexval($exon);
			    }
#			    print LEx $scoreag[$i]," ",$length-3," ",$hexv," 2\n";
			    print LExAll $name," ",$scoreag[$i]," ",$length-3," ",$hexv," ",$ag[$i]+2+$length-3," 2\n";
			}
		    }
		
		    $last+=3;
		    my $temp=substr($process,3);
		    $process=$temp;
		}
	    }
	    else {
# 		$length=$ex[$#ex]-$ex[$#ex-1]+1-3;
# 		if($length<6) { $hexv=0;}
# 		else {
# 		    $exon=substr($seq,$lastag+1,$length);
# 		    $hexv=hexval($exon);
# 		}
# 		print LEx $scoreag[$i]," ",$length," ",$hexv," 1\n";
	    }

	    # process internal exons
	    for(my $j=0;$j<=$igt;$j++) {
	    
		if($firstexon) {
		    # check here for first exons 
		    if($firstgt == $gt[$j]) {
# 			$length=$ex[2]-$ex[1]+1;
# 			if($length<6) { $hexv=0;}
# 			else {
# 			    $exon=substr($seq,$ex[1]-1,$length);
# 			    $hexv=hexval($exon);
# 			}
# 			print FEx $scoregt[$j]," ",$length," ",$hexv," 1\n";
		    }
		    else {
			$orf0=1;
			$orf1=1;
			$orf2=1;
			my $lastatg0=0;
			my $lastatg1=0;
			my $lastatg2=0;
			$process=substr($seq,0,$gt[$j]-1);
			$last=$gt[$j]-4;
			while($orf0||$orf1||$orf2) {
			    $length=$gt[$j]-2-$last+1;
			    if($length > $max_ex || $last-3<0) { last; }
			    if($orf0) {
				my $cod0=substr($process,-3);
				if(($cod0 eq "taa")||($cod0 eq "tag")||($cod0 eq "tga")) { 
				    $orf0=0;
				}
				if($cod0 eq "atg") {
				    $lastatg0=$last;
				}
				$last--;
				chop($process);
			    }
			    if($orf1) {
				my $cod1=substr($process,-3);
				if(($cod1 eq "taa")||($cod1 eq "tag")||($cod1 eq "tga")) { 
				    $orf1=0;
				}
				if($cod1 eq "atg") {
				    $lastatg1=$last;
				}
				$last--;
				chop($process);
			    }
			    if($orf2) {
				my $cod2=substr($process,-3);
				if(($cod2 eq "taa")||($cod2 eq "tag")||($cod2 eq "tga")) { 
				    $orf2=0;
				}
				if($cod2 eq "atg") {
				    $lastatg2=$last;
				}
				$last--;
				chop($process);
			    }
			}
			if($lastatg0) {
			    $length=$gt[$j]-2-$lastatg0+1;
			    if($length<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$lastatg0,$length);
				$hexv=hexval($exon);
			    }
#			    print FEx $scoregt[$j]," ",$length," ",$hexv," 2\n";
			    print FExAll $name," ",$scoregt[$j]," ",$length," ",$hexv," ",$lastatg0+1," 2\n";
			}
			if($lastatg1) {
			    $length=$gt[$j]-2-$lastatg1+1;
			    if($length<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$lastatg1,$length);
				$hexv=hexval($exon);
			    }
#			    print FEx $scoregt[$j]," ",$length," ",$hexv," 2\n";
			    print FExAll $name," ",$scoregt[$j]," ",$length," ",$hexv," ",$lastatg1+1," 2\n";
			}
			if($lastatg2) {
			    $length=$gt[$j]-2-$lastatg2+1;
			    if($length<6) { $hexv=0;}
			    else {
				$exon=substr($seq,$lastatg2,$length);
				$hexv=hexval($exon);
			    }
#			    print FEx $scoregt[$j]," ",$length," ",$hexv," 2\n";
			    print FExAll $name," ",$scoregt[$j]," ",$length," ",$hexv," ",$lastatg2+1," 2\n";
			}
		    }
		}

	    
		# now process the internal exons (only if it makes sense)
		if($gt[$j]>$ag[$i]+4) {
		    my $e1=$ag[$i]+2;
		    my $e2=$gt[$j]-1;
		    my $exonedges=" ".$e1." ".$e2." ";
		    $length = $gt[$j]-1-$ag[$i]-1;
		    if($length<=$max_ex) {
			if($length<6) { $hexv=0;}
			else {
			    $exon=substr($seq,$ag[$i]+1,$length);
			    $hexv=hexval($exon);
			}
			if(index($exons,$exonedges)>-1) { # true exon
			    print IEx $scoreag[$i]," ",$scoregt[$j]," ",$length," ",$hexv," 1\n";
			}
			else {
			    if(oneorf($exon) && $length>=6) {
				print IEx $scoreag[$i]," ",$scoregt[$j]," ",$length," ",$hexv," 2\n";
			    }
			}
		    }
		}
	    
		# process introns too
# 		if($gt[$j]+3<$ag[$i]) {
# 		    my $e1=$gt[$j]-1;
# 		    my $e2=$ag[$i]+2;
# 		    my $exonedges=" ".$e1."  ".$e2." ";
# 		    $length=$ag[$i]+1-$e1;
# 		    if($length<=$max_in) {
# 			if($length<6) { $hexv=0;}
# 			else {
# 			    $intron=substr($seq,$gt[$j]-1,$length);
# 			    $hexv=hexval($intron);
# 			}
# 			if(index($exons,$exonedges)>-1) { # true intron
# 			    print Intr $scoregt[$j]," ",$scoreag[$i]," ",$length," ",$hexv," 1\n";
# 			}
# 			else {
# 			    print Intr $scoregt[$j]," ",$scoreag[$i]," ",$length," ",$hexv," 2\n";
# 			}
# 		    }
# 		}
	    }
	
	    $firstexon=0;
	}
    }
}


# this procedure takes a feature file and transforms it into a balanced
# one (i.e. reduces the no. of features in class 2 by randomly 
# dropping some of them)
sub balance {

    my ($f,$o)=@_;

    my ($n1,$n2);

    open(F,$f) or die "ERROR 7: Couldn't open the features file $f\n";
    while(<F>) {
	chomp;
	my @a=split;
	if($a[$#a]==1) { $n1++;}
	if($a[$#a]==2) { $n2++;}
    }
    close(F);

    if(!$n1) { 
	die "ERROR 146: Feature file $f has not any true examples! Impossible to train the gene finder.\n";
    }

    my $droprate=int($n2/$n1);

    open(F,$f) or die "ERROR 8: Couldn't open the features file $f\n";
    open(O,">$o") or die "ERROR 9: Couldn't create the features file $o for the DT training\n";

    while(<F>){
	chomp;
	my @a=split;
	my $roll=0;
	if($a[$#a]==2) {
	    $roll=int(rand $droprate);
	}
	if(!$roll || $n2<$n1) {
	    print O $_,"\n";
	}
    }

    close(F);
    O->autoflush(1);
    close(O);

}


sub createtree_old {

    my ($class, $treename, $scriptdir)=@_;
    
    my (@max);

    my $dtname=$class.".dt";

    for(my $i=1;$i<=$no_of_trees;$i++) {
    	$max[$i]=0;
    }

    for(my $i=1;$i<300;$i++) {
    
	my $status=system("$scriptdir/mktree -t $class -p0.15 -s$i > out");
	die "ERROR 10: Could not execute program $scriptdir/mktree!\n" unless $status==0;

	open(F,"out") or die "ERROR 11: Could not open temporary decision tree output file!\n";

	while(<F>) {
	
	    chomp;
	    my @arg=split(/\s+/,$_);

	    if($arg[0] eq "acc.") {
	
		my $leaves=$arg[8];

		#print "$i: no_of_leaves=$leaves\n";

		if($leaves<=$no_leaves) {
	    
		    my $val=$arg[5];

		    my $equal=0;
		    my $j=$no_of_trees;
		    while($j>0) {
			if($val>$max[$j]) { $j=0; }
			elsif($val==$max[$j]) { 
			    $equal=1;
			    $j=0;
			}
			else {$j--;}
			
		    }
		    
		    $j=1;
		    if(!$equal) {

			my $minacc=101;
			my $minval=0;
			for($j=1;$j<=$no_of_trees;$j++) {
			    if($val>$max[$j]) {
				if($max[$j]<$minacc) {
				    $minacc=$max[$j];
				    $minval=$j;
				}
			    }
			}
			if($minval) {
			    print "$i: ";
			    system("cat out");
			    
			    $status=system("cp $dtname trees/$treename."."dt.$minval");
			    die "ERROR 13: Could not copy $dtname to trees/." unless $status==0;
			    $max[$minval]=$val;
			}

		    }
		}
	    }
	}

	close(F);
    }

}

sub createtree {

    my ($class, $treename, $scriptdir)=@_;
    
    my (@max,@leaf);

    my $dtname=$class.".dt";

    for(my $i=1;$i<=$no_of_trees;$i++) {
    	$max[$i]=0;
    }

    for(my $i=1;$i<300;$i++) {
    
	my $status=system("$scriptdir/mktree -a -t $class -p0.15 -s$i > out");
	die "ERROR 10: Could not execute program $scriptdir/mktree!\n" unless $status==0;

	open(F,"out") or die "ERROR 11: Could not open temporary decision tree output file!\n";

	while(<F>) {
	
	    chomp;
	    my @arg=split(/\s+/,$_);

	    if($arg[0] eq "acc.") {
	
		my $leaves=$arg[8];

		#print "$i: no_of_leaves=$leaves\n";

		#if($leaves<=$no_leaves) {
	    
		    my $val=$arg[5];

		    my $equal=0;
		    my $minval=0;
		    my $j=$no_of_trees;
		    while($j>0) {
			if(!($leaf[$j])) { $minval=$j; $j--; }
			elsif($leaves<$leaf[$j]) { $minval=$j; $j--;}
			elsif($leaves==$leaf[$j]) { 
			    if($val<=$max[$j]) {
				$equal=1;
				$j=0;
			    }
			    else { $minval=$j;$j--;}
			}
			else {$j--;}
			
		    }
		    
		    $j=1;
		    if(!$equal && $minval) {

			my $minacc=1;
			   $minval=0;
			for($j=1;$j<=$no_of_trees;$j++) {
			    if(!($leaf[$j])) {$minval=$j; last;}
			    if($leaves==$leaf[$j]) { $minval=$j; last;}
			    if($leaves<$leaf[$j]) {
				if($leaf[$j]>$minacc) {
				    $minacc=$leaf[$j];
				    $minval=$j;
				}
			    }
			}
			if($minval) {
			    print "$i: ";
			    system("cat out");
			    #print "$i: copy to trees/$treename.dt.$minval\n";
			    $status=system("cp $dtname trees/$treename."."dt.$minval");
			    die "ERROR 13: Could not copy $dtname to trees/." unless $status==0;
			    $max[$minval]=$val;
			    $leaf[$minval]=$leaves;

			}

		    }
		#}
	    }
	}

	close(F);
    }
    
    my $k=0;
    my $maxacc=0;
    my $maxval=0;
    for(my $j=1;$j<=$no_of_trees;$j++) {
	if(!$leaf[$j]) { if(!$k) { $k=$j; } }
	else {
	    if($max[$j]>$maxacc) {
		$maxacc=$max[$j];
		$maxval=$j;
	    }
	}
    }
    if(!$maxval) {
	die "ERROR 147: Not enough examples to create the $treename trees!\n";
    }
    else {
	if($k) {
	    for(my $j=$k;$j<=$no_of_trees;$j++) {
		my $status=system("cp trees/$treename."."dt.$maxval trees/$treename."."dt.$j");
		die "ERROR 12: Could not copy trees/$treename.dt.$maxval" unless $status==0; 
	    }
	}
    }
}


sub printtreenames {
    my ($output,$name)=@_;

    open(F,">>$output") or die "ERROR 14: Could not open $output file for writing!\n";

    for(my $i=1;$i<11;$i++) {
	print F "$name.dt.$i\n";
    }

    close(F);
}

# this procedure extracts exons and introns and
# prepares them for IMM scoring
# usage dtorfs.pl seqs exons.dat sites.scores
sub dtorfs {
    my ($seqs,$exons,$sites,$dtorfsfile) = @_;
    my %ex;

    open(F,$exons);
    
    while(<F>) {
	chomp;
	if($_) {
	    my @a=split;
	    $ex{$a[0]}.=$a[1]." ".$a[2]." ";
	}
    }
    close(F);

    my ($gene,%dgt,%dag,%dscoregt,%dscoreag);

    open(F,$sites);
    while(<F>) {
	chomp;
	my @a=split;
	if($a[0] eq "Gene") { $gene=$a[1];}
	elsif($a[0] eq "gt:") {
	    my $pos=$a[1]-1;
	    $dgt{$gene}.=$pos." ";
	    $dscoregt{$gene}.=$a[2]." ";
	}
	elsif($a[0] eq "ag:") {
	    my $pos=$a[1]+2;
	    $dag{$gene}.=$pos." ";
	    $dscoreag{$gene}.=$a[2]." ";
	}
    }
    close(F);

    open(F,$seqs);
    open(D,">$dtorfsfile");
    while(<F>) {
	chomp;
	my ($name,$seq)=/^(\S+)\s*(\S+)$/;
    
	#print "$name ",$ex{$name},"\n";

	my @e=split(/\s+/,$ex{$name});
	@gts=split(/\s+/,$dgt{$name});
	@ags=split(/\s+/,$dag{$name});
	@sgt=split(/\s+/,$dscoregt{$name});
	@sag=split(/\s+/,$dscoreag{$name});

	my ($len,$seqex,$seqin);

	my ($score1,$score2);

	for(my $i=0;$i<=$#e;$i+=2){
	    
	    $len=$e[$i+1]-$e[$i]+1;
	    $seqex=substr($seq,$e[$i]-1,$len);
	    if($i+2>$#e) { chop($seqex);chop($seqex);chop($seqex); }

	    #print "$name ",$ex{$name}," $i $seqex\n";exit;

	    if($#e==1 && $seqex) {
		print D "S $name ",$i/2+1," $seqex\n";
	    }
	    else{
		if($i==0 && $seqex) {
		    $score2=search_score($e[$i+1],2);
		    print D "F $name ",$i/2+1," $seqex $score2\n";
		}
		else {
		    $score1=search_score($e[$i],1);
		    $seqin=substr($seq,$e[$i-1],$e[$i]-1-$e[$i-1]);
		    if($seqin) { print D "I $name ",$i/2+1," $seqin $score2 $score1\n";}
		    if($i+1==$#e && $seqex) {
			print D "L $name ",$i/2+1," $seqex $score1\n";
		    }
		    else {
			if($seqex) {
			    $score2=search_score($e[$i+1],2);
			    print D "E $name ",$i/2+1," $seqex $score1 $score2\n";
			}
		    }
		}
	    }
	}
    }
    close(F);
    close(D);
}

sub search_score {
    my ($pos,$type)=@_;
    
    my $score;
    
    if($type==1) {
	my $i=0;
	$score=-99;
	while($i<=$#ags) {
	    if($pos==$ags[$i]) { return($sag[$i]);}
	    $i++;
	}
    }
    
    if($type==2) {
	my $i=0;
	$score=-99;
	while($i<=$#gts) {
	    if($pos==$gts[$i]) { return($sgt[$i]);}
	    $i++;
	}
    }

    return($score);
}


# this procedure extracts the features file for the upstream
# and downstream regions of the coding gene
sub features_updown {
    my ($seqs,$satg,$sstop,$umodel,$dmodel,$cmodel,$scriptdir,$freq_a,$freq_c,$freq_g,$freq_t)=@_;

    open(F,$satg) ||  die "ERROR 169: Could not open file $satg for reading: $!\n"; # sites.atg
    <F>;
    my $true=1;
    my %truestart;
    my %falsestart;
    while(<F>) {
	chomp;
	if($_) {
	    if($_ eq "Scores for false starts") { $true=0; }
	    else {
		my @a=split;
		if($true) {
		    $truestart{$a[0]}=$a[1]." ".$a[2];
		}
		else { 
		    $falsestart{$a[0]}.=$a[1]." ".$a[2]." ";
		}
	    }
	}
    }
    close(F);


    open(F,$sstop) ||  die "ERROR 170: Could not open file $sstop for reading: $!\n"; # sites.stop
    <F>;
    $true=1;
    my %truestop;
    my %falsestop;
    while(<F>) {
	chomp;
	if($_) {
	    if($_ eq "Scores for false stops") { $true=0; }
	    else {
		my @a=split;
		if($true) {
		    $truestop{$a[0]}=$a[1]." ".$a[2];
		}
		else { 
		    $falsestop{$a[0]}.=$a[1]." ".$a[2]." ";
		}
	    }
	}
    }
    close(F);

    my  $upfile="upstream.features";
    my  $downfile="downstream.features";

    open(U,">$upfile");
    open(D,">$downfile");

    my $sum_all=$freq_a+$freq_c+$freq_g+$freq_t;
    $freq_a/=$sum_all;
    $freq_c/=$sum_all;
    $freq_g/=$sum_all;
    $freq_t/=$sum_all;

    open(F,$seqs) || die "ERROR 171: Could not open file $seqs for reading: $!\n";    # seqs
    while(<F>) {
	chomp;
	my ($name,$seq) = /^(\S+)\s+(\S+)$/;
	
	my $len=length($seq);

	if($truestart{$name} || $falsestart{$name}) {
	    my @t;
	    my @f;

	    @t=split(/\s+/,$truestart{$name});
	    @f=split(/\s+/,$falsestart{$name});

	    if($truestart{$name} && $t[0] >= 100) { 
		my $secv=substr($seq,$t[0]-100,99); 

		#print "$scriptdir/updomeasure $secv 9 90 $umodel $dmodel $cmodel 1 $freq_a $freq_c $freq_g $freq_t ; t[0]=",$t[0]," ",$truestart{$name};

		my $output=`$scriptdir/updomeasure $secv 9 90 $umodel $dmodel $cmodel 1 $freq_a $freq_c $freq_g $freq_t`; # $output contains the two values 
		if($?) { die printerr("ERROR 172: $scriptdir/updomeasure exited funny: $?");}
		print U $t[1]," $output 1\n";
	    }
	    if($falsestart{$name} && $f[0] >= 100) {
		my $secv=substr($seq,$f[0]-100,99); 

		#print "$scriptdir/updomeasure $secv 9 90 $umodel $dmodel $cmodel 1 $freq_a $freq_c $freq_g $freq_t ; f[0]=",$f[0]," ",$falsestart{$name};

		my $output=`$scriptdir/updomeasure $secv 9 90 $umodel $dmodel $cmodel 1 $freq_a $freq_c $freq_g $freq_t`; # $output contains the two values 
		if($?) { die printerr("ERROR 173: $scriptdir/updomeasure exited funny: $?"); }
		print U $f[1]," $output 2\n";
	    }
	}
	
	#print "$name: ",$truestop{$name}," ",$falsestop{$name},"\n";

	if($truestop{$name} || $falsestop{$name}) {
	    my @t;
	    my @f;
	    @t=split(/\s+/,$truestop{$name});
	    @f=split(/\s+/,$falsestop{$name});

	    if($truestop{$name} && $len-$t[0] >= 90) { 
		my $secv=substr($seq,$t[0]-11,102); 

		#print "$scriptdir/updomeasure $secv 11 90 $umodel $dmodel $cmodel -1 $freq_a $freq_c $freq_g $freq_t ; t[0]=",$t[0]," ",$truestop{$name};

		my $output=`$scriptdir/updomeasure $secv 11 90 $umodel $dmodel $cmodel -1 $freq_a $freq_c $freq_g $freq_t`; # $output contains the two values 
		if($?) { die printerr("ERROR 174: $scriptdir/updomeasure exited funny: $?");}
		print D $t[1]," $output 1\n";
	    }
	    if($falsestop{$name} && $len-$f[0] >= 90) {
		my $secv=substr($seq,$f[0]-11,102);

		#print "$scriptdir/updomeasure $secv 11 90 $umodel $dmodel $cmodel -1 $freq_a $freq_c $freq_g $freq_t ; f[0]=",$f[0]," ",$falsestop{$name};
 
		my $output=`$scriptdir/updomeasure $secv 11 90 $umodel $dmodel $cmodel -1 $freq_a $freq_c $freq_g $freq_t`; # $output contains the two values 
		if($?) { die printerr("ERROR 175: $scriptdir/updomeasure exited funny: $?");}
		print D $f[1]," $output 2\n";
	    }
	}
    }
 
    close(F);
    close(U);
    close(D);
}
    
# oneorf is used by feature_extract 
#sub oneorf {
#    my $prelex=$_[0];
#    my $ret=3;
#    my @frame;
#    my @test=split(/tga|taa|tag/,$prelex);
#    my $length=0;
#    my $i=0;
    
#    $frame[0]=1;
#    $frame[1]=1;
#    $frame[2]=1;

#    while($ret && $i<=$#test) {
#	$length+=length($test[$i])+3;
#	my $val=$length % 3;
#	$frame[$val]=0;
#	$ret=$frame[0]+$frame[1]+$frame[2];
#	$i++;
#    }

#    return($ret);
#}

sub DTtraintrue {

    # input exons.dat seqs sites.scores sites.atg sites.stop

    my ($exfile,$seqs,$splicesc,$startsc,$stopsc,$fexfile,$lexfile,$sexfile)=@_;

    # read exons
    open(F,$exfile);
    my %exons;
    while(<F>) {
        chomp;
        if($_) {
	    my ($name,$beg,$end)=split;
	    $exons{$name}.=$beg." ".$end." ";
        }
    }
    close(F);


    open(F,$splicesc) || die "ERROR 176: Could not open file $splicesc for reading: $!\n"; # sites.scores

    my ($iex,$name,@e);

    my (%gtid,%agid);

    while(<F>) {
        chomp;
        my @a=split;
        my $id;

        if($a[0] eq "Gene") { 
	    $name=$a[1];
	    @e=split(/\s+/,$exons{$name});
	    $iex=1;
        }
        else {
	    if($iex<=$#e && $a[0] eq "gt:" && $iex%2==1 && $a[1]==$e[$iex]+1) {
	        $id=$name." ".$e[$iex];
	        $gtid{$id}=$a[2];
	        $iex++;
	    }   
	    elsif($iex<=$#e && $a[0] eq "ag:" && $iex%2==0 && $a[1]==$e[$iex]-2) {
	        $id=$name." ".$e[$iex];
	        $agid{$id}=$a[2];
	        $iex++;
	    }
	    elsif($iex<=$#e && $a[1]>$e[$iex]+1) { $iex++;}
        }
    }

    close(F);

    open(F,$startsc) || die "ERROR 177: Could not open file $startsc for reading: $!\n"; # sites.atg

    my $true;
    my %starttrue;
    while(<F>) {
        chomp;
        my @a=split;

        if($a[2] eq "true") { $true=1;}
        elsif($a[2] eq "false") { $true=0;}
        elsif($true) {
	    $starttrue{$a[0]}=$a[2];
        }
    }
    close(F);

    open(F,$stopsc) || die "ERROR 178: Could not open file $stopsc for reading: $!\n"; # sites.stop

    my %stoptrue;
    while(<F>) {
        chomp;
        my @a=split;

        if($a[2] eq "true") { $true=1;}
        elsif($a[2] eq "false") { $true=0;}
        elsif($true) {
	    $stoptrue{$a[0]}=$a[2];
        }
    }
    close(F);

    # print features
    if(!$fexfile) { $fexfile="tfirst.train";}
    if(!$lexfile) { $lexfile="tlast.train";}
    if(!$sexfile) { $sexfile="tsgle.train";}

    open(FEx,">$fexfile");
#    open(IEx,">$iexfile");
    open(LEx,">$lexfile");
    open(SEx,">$sexfile");

    open(F,$seqs);

    while(<F>) {
        chomp;
        my ($name,$seq)=split;

        @e=split(/\s+/,$exons{$name});
    
        if($#e==1) { # single exon
	    my $len=$e[1]-$e[0]+1;
	    my $truex=substr($seq,$e[0]-1,$len);
	    my $hexv=hexval($truex,$len);
	    # print start_site_score stop_site_score length hexamerval 2-iftrue
	    if($starttrue{$name} && $stoptrue{$name}) { 
                print SEx $starttrue{$name}," ",$stoptrue{$name}," ",$len," ",$hexv," 1\n";
            }
        }
        else { # all other exons
	    for(my $i=0;$i<$#e;$i+=2) {
	        my $len=$e[$i+1]-$e[$i]+1;
	        my $truex=substr($seq,$e[$i]-1,$len);
	        my $hexv=hexval($truex,$len);

	        if($i==0) {
		    my $id2=$name." ".$e[$i+1];
		    if($starttrue{$name} && $gtid{$id2} ) { 
                        print FEx $starttrue{$name}," ",$gtid{$id2}," ",$len," ",$hexv," 1\n";
                    }
	        }
	        elsif($i==$#e-1) {
		    my $id1=$name." ".$e[$i];
		    if($agid{$id1} && $stoptrue{$name}) { 
                        print LEx $agid{$id1}," ",$stoptrue{$name}," ",$len," ",$hexv," 1\n";
                    }
	        }
	        else {
		    my $id1=$name." ".$e[$i];
		    my $id2=$name." ".$e[$i+1];
		    if($agid{$id1} && $gtid{$id2}) { 
                        print IEx $agid{$id1}," ",$gtid{$id2}," ",$len," ",$hexv," 1\n";
                    }
	        }
	    }
        }
    }
    close(F);
    close(FEx);
    close(IEx);
    close(LEx);
    close(SEx);
}

sub allinfo2first {
    my ($sitesfalse,$allinfofeat,$outpf)=@_;

    # first "sites.false";

    open(F,$sitesfalse);
    my %atg;
    while(<F>) {
        chomp;
        my @a=split;
        if($a[1] eq "atg") {
	    my $name=$a[0]."_".$a[2];
	    $atg{$name}=$a[3];
        }
    }
    close(F);

    # second "firstexons.features.allinfo";

    open(F,$allinfofeat);
    open(O,">$outpf");
    while(<F>) {
        chomp;
        my @a=split;
        my $name=$a[0]."_".$a[4];
        if($atg{$name}) {
	    print O $atg{$name}," ",$a[1]," ",$a[2]," ",$a[3]," 2\n";
        }
    }
    close(F);
}

sub allinfo2last {
    my ($sitesfalse,$allinfofeat,$outpf)=@_;

    # first "sites.false";

    open(F,$sitesfalse);
    my %stop;
    while(<F>) {
        chomp;
        my @a=split;
        if($a[1] eq "stop") {
	    my $name=$a[0]."_".$a[2];
	    $stop{$name}=$a[3];
        }
    }
    close(F);

    # second "lastexons.features.allinfo";

    open(F,$allinfofeat);
    open(O,">$outpf");
    while(<F>) {
        chomp;
        my @a=split;
        my $last=$a[4];
        my $name=$a[0]."_".$last;
        if($stop{$name}) {
	    print O $a[1]," ",$stop{$name}," ",$a[2]+3," ",$a[3]," 2\n";
        }
    }
    close(F);
}   

sub allinfo2sgle {
    my ($sitesfalse,$allinfofeat,$outpf)=@_;

    # first "sites.false";

    open(F,$sitesfalse);
    my (%atg,%stop);
    while(<F>) {
        chomp;
        my @a=split;
        if($a[1] eq "atg") {
	    my $name=$a[0]."_".$a[2];
	    $atg{$name}=$a[3];
        }
        if($a[1] eq "stop") {
	    my $name=$a[0]."_".$a[2];
	    $stop{$name}=$a[3];
        }
    }
    close(F);

    # second "snglexons.features.allinfo";

    open(F,$allinfofeat);
    open(O,">$outpf");
    while(<F>) {
        chomp;
        my @a=split;
        my $name1=$a[0]."_".$a[3];
        my $last=$a[1]+$a[3];
        my $name2=$a[0]."_".$last;
        if($atg{$name1}&&$stop{$name2}) {
	    print O $atg{$name1}," ",$stop{$name2}," ",$a[1]," ",$a[2]," 2\n";
        }
    }
    close(F);
}
