#!/usr/local/bin/perl

#Copyright (c) 2003 by Mihaela Pertea.

# the scripts computes statistics for the cfg file

use strict;
use FileHandle;

my $initin0=0.333;
my $initin1=0.333;
my $initin2=0.333;

return 1;

sub printstatistics {

    my ($cds,$outpf,$mean5utr,$mean3utr,$meaninterg)=@_; 

    my($name,$beg,$end,$lastbeg,$lastend);
    my $noex;
    my $nogene=0;
    my $nosingle=0;
    my ($phase,$lastphase,$length);
    my $intronno=0;
    my $intronlen=0;
    my (%start,%stop,%seen);
    my (@pinitin,@pinin);
    my $pintrex=0;
    my $pterm=0;
    my $exlength=0;
    my $allexon=0;

    open(F,$cds);

    while(<F>) {
        chomp;
        if($_) {
	    my ($name,$beg,$end)=split;
	    $exlength+=$end-$beg+1;
	    $allexon++;
	    if(!$seen{$name}) { 
	        $start{$name}=$beg;
	        $seen{$name}++;
	        $noex=1;
	        $length=$end-$beg+1;
	        $nogene++;
	    }   
	    else {
		$phase=$length%3;
		if($noex==1) {
		    $pinitin[$phase]++;
		}
		else {
		    $pinin[$lastphase][$phase]++;
		}
		$intronlen+=$beg-1-$lastend;
		$intronno++;
		$noex++;
		$length+=$end-$beg+1;
		$lastphase=$phase;
		$pintrex++;
	    }
	    $lastbeg=$beg;
	    $lastend=$end;
	}
	else {
	    $stop{$name}=$lastend;
	    if($noex>1) {
		$pintrex--;
		$pterm++;
	    }
	    else {
		$nosingle++;
	    }
	}
    }
    close(F);

    ### this is very approximative; 
    my $exavglen=$exlength/$nogene;
    my $inavglen;
    if($intronno) { $inavglen=$intronlen/$intronno; }
    else { $inavglen=0;}
    my $inavgpergene=$inavglen*($allexon/$nogene-1);
    my $allavglen=$mean5utr+$mean3utr+$meaninterg+$exavglen+$inavgpergene;
    if($mean5utr) { printf $outpf "Init5'UTR %.6f\n",$mean5utr/$allavglen };
    if($mean3utr) { printf $outpf "Init3'UTR %.6f\n",$mean3utr/$allavglen };
    if($meaninterg) { 
	printf $outpf "InitInterg %.6f\n",$meaninterg/$allavglen;
	printf $outpf "InitIntron %.6f\n",$inavgpergene/$allavglen;
    }
    else {
	$meaninterg=250;
	printf $outpf "InitInterg 0.95\n";
	printf $outpf "InitIntron 0.05\n";
    }
    ###

    printf $outpf "InitI0 %.6f\n",$initin0;
    printf $outpf "InitI1 %.6f\n",$initin1;
    printf $outpf "InitI2 %.6f\n",$initin2;
    printf $outpf "MeanIntergen %.1f\n",$meaninterg;
    printf $outpf "Mean5'UTR %.1f\n",$mean5utr if($mean5utr);
    printf $outpf "MeanIntron %.1f\n",$inavglen;
    printf $outpf "Mean3'UTR %.1f\n",$mean3utr if($mean3utr);
    printf $outpf "PSngl %.6f\n",$nosingle/$nogene;
    if($pterm+$pintrex) { printf $outpf  "PTerm %.6f\n",$pterm/($pterm+$pintrex); }
    else { printf $outpf  "PTerm %.6f\n",0;}
    my $sum=$pinitin[0]+$pinitin[1]+$pinitin[2];
    for(my $i=0;$i<3;$i++) {
	if($sum) { printf $outpf "PInitIn%d %.3f\n",$i,$pinitin[$i]/$sum; }
	else { printf $outpf "PInitIn%d %.3f\n",$i,0;}
    }
    for(my $i=0;$i<3;$i++) {
	$sum=$pinin[$i][0]+$pinin[$i][1]+$pinin[$i][2];
	for(my $j=0;$j<3;$j++) {
	    if($sum) { printf $outpf "PIn%dIn%d %.3f\n",$i,$j,$pinin[$i][$j]/$sum;}
	    else { printf $outpf "PIn%dIn%d %.3f\n",$i,$j,0;}
	}
    }
}

sub lendistr {
    my ($exons,$outf,$nval,$smoothwindow,$erfappdir)=@_;

    open(F,$exons);
    my $lastname="";
    my $no_ex=0;
    my (@ex);
    my $lastlen=-1;

    while(<F>) {
	chomp;
	if($_) {
	    my ($name,$beg,$end)=split;
	    if($name ne $lastname) { # new label
		if($no_ex>1) { # last exon
		    $ex[2][$lastlen]++;
		}
		elsif($no_ex==1) { # single exon
		    $ex[3][$lastlen]++;
		}
		$no_ex=1;
		$lastname=$name;
	    }
	    else { # same gene
		if($no_ex==1) { # first exon
		    $ex[0][$lastlen]++;
		}
		else { # internal exon
		    $ex[1][$lastlen]++;
		}
		$no_ex++;
	    }
	    $lastlen=int(($end-$beg+1)/3);
	    $lastlen=$nval-1 if($lastlen>=$nval);
	}
    }

    close(F);
    
    if($no_ex>1) { # last exon
	$ex[2][$lastlen]++;
    }
    elsif($no_ex==1) { # single exon
	$ex[3][$lastlen]++;
    }

    open(F,">$outf");

    for(my $i=0;$i<4;$i++) {
	if(scalar($ex[$i])) { smooth($nval,$ex[$i],$smoothwindow,$erfappdir);}
	if($i==0) { 
	    print F "Initial $nval\n";
	}
	elsif($i==1) {
	    print F "Internal $nval\n";
	}
	elsif($i==2) {
	    print F "Terminal $nval\n";
	}
	elsif($i==3) {
	    print F "Single $nval\n";
	}

	for(my $j=0;$j<$nval;$j+=5) {
	    printf F "%.6f        %.6f        %.6f        %.6f        %.6f\n",$ex[$i][$j],$ex[$i][$j+1],$ex[$i][$j+2],$ex[$i][$j+3],$ex[$i][$j+4];
	}
    }

    close(F);
}

sub gfunc{
    my ($j,$mean,$sigma,$erfappdir)=@_;

    my $val=0;

    my $j1=($j+0.5-$mean)/($sigma*sqrt(2));
    my $x1=`$erfappdir/erfapp $j1`;

    my $j2=($j-0.5-$mean)/($sigma*sqrt(2));
    my $x2=`$erfappdir/erfapp $j2`;
    
    $val=($x1-$x2)/2;
}


sub smooth {
    my ($n,$lengths,$smoothwindow,$erfappdir)=@_;
    my $eps=0.000001;
    my $N=0;
    my @sum;

    for(my $k=0;$k<$n;$k++) { $N+=$$lengths[$k];}

    for(my $x=0;$x<$n;$x++) {
	
	$sum[$x]=0;

	my $k = $x-$smoothwindow > 0 ? $x-$smoothwindow : 0;
	while($k<$n && $k<=$x+$smoothwindow) {
	    if($$lengths[$k]>0) {
		my $sigmak=sqrt(2*($k+1)/$$lengths[$k]);
		$sum[$x]+=$$lengths[$k]*gfunc($x+1,$k+1,$sigmak,$erfappdir);
	    }
	    $k++;
	}

	$sum[$x]/=$N;
    }

    my $allsum=0;

    for(my $x=0;$x<$n;$x++) {
	$allsum+=$sum[$x];
    }

    $allsum/=(1-$n*$eps);

    for(my $k=0;$k<$n;$k++) {
	$$lengths[$k]=$eps+$sum[$k]/$allsum;
    }

}



