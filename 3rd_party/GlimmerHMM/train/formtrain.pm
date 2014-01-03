#!/usr/local/bin/perl
#Copyright (c) 2003 by Mihaela Pertea.

use strict;
use FileHandle;

return 1;

sub formacc {
    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 101: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 102: Couldn't open $ex for reading: $!\n";

    my $ind=0;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=0;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start,$end}++;
	    }
	    $ind=1;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 103: Couldn't open $seq for reading: $!\n";

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;
	my $k;

	foreach $k (keys %isexon)
	{
	    my ($acc,$start,$end)=split(/$;/,$k);
	    if($acc eq $name)
	    {
		my $secv=substr($s,$start-71,80);
		print O "Acceptor $secv $name\n";
	    }
	}
    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub formfacc {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 104: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 105: Couldn't open $ex for reading: $!\n";

    my $ind=0;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=0;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start}++;
	    }
	    $ind=1;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 106: Couldn't open $seq for reading: $!\n";

    my (@poz,$val);

    while(<F>)
    {
	chomp;
	my ($name,$s)=split;
	my $l=length($s);

	my $count=0;
	for(my $i=70;$i<$l-12;$i++)
	{
	    if(substr($s,$i,2) eq "ag")
	    {
		if(!$isexon{$name,$i+3}) { $poz[$count++]=$i+3;}
	    }
	}

	for(my $i=0;$i<60;$i++)
	{
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count))
	    {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }

	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-71,80);
	    print O "FA $secv $name\n";
	    $val=0;
	}

    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub formdon {
    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 107: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 108: Couldn't open $ex for reading: $!\n";
    
    my (%isexon,$name1,$start1,$end1);

    while(<E>)
    {
	chomp;

	my ($name,$start,$end)=split;

	if($_ ne "")
	{
	    $isexon{$name1,$start1,$end1}++;
	}

	$name1=$name;
	$start1=$start;
	$end1=$end;
    }

    close(E);

    open(F,$seq) or die "ERROR 109: Couldn't open $seq for reading: $!\n";

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;
	
	my $k;
	foreach $k (keys %isexon)
	{
	    my ($acc,$start,$end)=split(/$;/,$k);
	    if($acc eq $name)
	    {
		my $secv=substr($s,$end-10,80);
		print O "Donor $secv $name\n";
	    }
	}
    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub formfdon {
    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 110: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 111: Couldn't open $ex for reading: $!\n";

    my (%isexon,$name1,$start1,$end1);
    
    while(<E>)
    {
	chomp;

	my ($name,$start,$end)=split;

	if($_ ne "")
	{
	    $isexon{$name1,$end1}++;
	}

	$name1=$name;
	$start1=$start;
	$end1=$end;
    }

    close(E);

    open(F,$seq) or die "ERROR 112: Couldn't open $seq for reading: $!\n";

    my (@poz,$val);

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;

	my $l=length($s);

	my $count=0;
	for(my $i=11;$i<$l-70;$i++)
	{
	    if(substr($s,$i,2) eq "gt")
	    {
		if(!$isexon{$name,$i}) { $poz[$count++]=$i;}
	    }
	}

	for(my $i=0;$i<60;$i++)
	{
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count))
	    {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }
	    
	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-10,80);
	    print O "FD $secv $name\n";
	    $val=0;
	}
	
    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub clean {
    my ($f,$o,$output)=@_;
    open(O,">$output") or die "ERROR 113: Couldn't open $output for writing: $!\n";

# $o is the option what to clean;

    open(F,$f) or die "ERROR 114: Couldn't open $f for reading: $!\n";

    my $atgf="atg.errors";
    my $accf="acc.errors";
    my $donf="don.errors";
    my $stopf="stop.errors";

    while(<F>){
	chomp;
	my @a=split;
	my $line=$_;
	
	my $l=length($a[1]);

	if($l==19) {
	    if($o eq "atg") {
		my $sir=substr($a[1],12,3);
		if(lc($sir) eq "atg") { print O $_,"\n";} 
		else { open(S,">> $atgf");print S "Error: $line\n";close(S);}
	    }
	    if($o eq "stop") {
		my $sir=substr($a[1],4,3);
		if((lc($sir) eq "taa")||(lc($sir) eq "tga")||(lc($sir) eq "tag")) { print O $_,"\n";} 
		else { open(S,">> $stopf");print S "Error: $line\n";close(S);}
	    }
	}

	if($l==80) {

	    if($o eq "don") {    
		my $sir=substr($a[1],10,2);
		
		if(lc($sir) eq "gt") { print O $_,"\n";} 
		else { open(D,">> $donf"); print D "Error: $line\n";close(D);}
	    }

	    if($o eq "acc") {
		my $sir=substr($a[1],68,2);
	    
		if(lc($sir) eq "ag") { print O $_,"\n";}
		else { open(A,">> $accf"); print A "Error: $line\n"; close(A);} 
	    }
	    if($o eq "len") { print O $_,"\n";}
	}
	if($l==40) {
	    if($o eq "gtag") {
		my $sir=substr($a[1],18,4);

		if(lc($sir) eq "aggt") { print O $_,"\n";} 
		else { print STDERR "Error sir=$sir: ",$_,"\n";}
	    }
	}
	
    }

    close (F);
    O->autoflush(1);
    close(O);
}

sub formacc162 {
    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 115: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 116: Couldn't open $ex for reading: $!\n";

    my $ind=0;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=0;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start,$end}++;
	    }
	    $ind=1;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 117: Couldn't open $seq for reading: $!\n";

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;

	my $k;
	foreach $k (keys %isexon)
	{
	    my ($acc,$start,$end)=split(/$;/,$k);
	    if($acc eq $name)
	    {
		my $secv=substr($s,$start-83,162);
#	    print "Acceptor $secv $name\n"; // initial variant
		printf O "A%ld $secv $name\n",$start-2;
	    }
	}
    }

    close(F);
    O->autoflush(1);
    close(O);
}

sub formfacc162 {
    
    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 118: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 119: Couldn't open $ex for reading: $!\n";

    my $ind=0;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=0;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start}++;
	    }
	    $ind=1;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 120: Couldn't open $seq for reading: $!\n";

    my (@poz,$val);

    while(<F>)
    {
	chomp;
	my ($name,$s)=split;
	my $l=length($s);

	my $count=0;
	for(my $i=100;$i<$l-100;$i++)
	{
	    if(substr($s,$i,2) eq "ag")
	    {
		if(!$isexon{$name,$i+3}) { $poz[$count++]=$i+3;}
	    }
	}

# to generate only a few false splice sites
	for(my $i=0;$i<30;$i++)
	{
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count))
	    {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }

	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-83,162);
#	$secv=substr($s,$val-$dist-3,2*$dist+2);
#	print "FA $secv $name\n";
	    printf O "FA%ld $secv $name\n",$val-2;
	    $val=0;
	}

    }

    close(F);
    O->autoflush(1);
    close(O);
}

sub formdon162 {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 121: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 122: Couldn't open $ex for reading: $!\n";

    my (%isexon,$name1,$start1,$end1);

    while(<E>)
    {
	chomp;

	my ($name,$start,$end)=split;
	
	if($_ ne "")
	{
	    $isexon{$name1,$start1,$end1}++;
	}

	$name1=$name;
	$start1=$start;
	$end1=$end;
    }

    close(E);

    open(F,$seq) or die "ERROR 123: Couldn't open $seq for reading: $!\n";
    
    while(<F>)
    {
	chomp;
	
	my ($name,$s)=split;

	my $k;
	foreach $k (keys %isexon)
	{
	    my ($acc,$start,$end)=split(/$;/,$k);
	    if($acc eq $name)
	    {
		my $secv=substr($s,$end-80,162);
#	    print "Donor $secv $name\n"; initial variant
		printf O "D%ld $secv $name\n",$end+1;
	    }
	}
    }
    
    close(F);
    O->autoflush(1);
    close(O);
}

sub formfdon162 {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 124: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 125: Couldn't open $ex for reading: $!\n";
    
    my (%isexon,$name1,$end1,$start1);

    while(<E>)
    {
	chomp;

	my ($name,$start,$end)=split;

	if($_ ne "")
	{
	    $isexon{$name1,$end1}++;
	}
	
	$name1=$name;
	$start1=$start;
	$end1=$end;
    }

    close(E);

    open(F,$seq) or die "ERROR 126: Couldn't open $seq for reading: $!\n";
    
    my (@poz,$val);

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;

	my $l=length($s);

	my $count=0;
	for(my $i=82;$i<$l-82;$i++)
	{
	    if(substr($s,$i,2) eq "gt")
	    {
		if(!$isexon{$name,$i}) { $poz[$count++]=$i;}
	    }
	}

	for(my $i=0;$i<30;$i++)
	{
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count))
	    {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }

	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-80,162);
	    #print "FD $secv $name\n";
	    printf O "FD%ld $secv $name\n",$val+1;
	    $val=0;
	}
	
    }
    
    close(F);
    O->autoflush(1);
    close(O);
}

sub clean162 {

    my ($f,$o,$output)=@_;
    open(O,">$output") or die "ERROR 127: Couldn't open $output for writing: $!\n";

# $o is the option what to clean;

    open(F,$f) or die "ERROR 128: Couldn't open $f for reading: $!\n";

    while(<F>){
	chomp;
	my @a=split;

	my $l=length($a[1]);

	if($l==162) {
	    
	    if($o eq "don") {    
		my $sir=substr($a[1],80,2);
	
		if($sir eq "gt") { print O $_,"\n";} 
		#else { print STDERR "Error: ",$_,"\n";}
	    }

	    if($o eq "acc") {
		my $sir=substr($a[1],80,2);
	    
		if($sir eq "ag") { print O $_,"\n";}
		#else { print STDERR "Error: ",$_,"\n";} 
	    }

	    if($o eq "gtag") {
		my $sir=substr($a[1],18,4);

		if($sir eq "aggt") { print O $_,"\n";} 
		#else { print STDERR "Error: ",$_,"\n";}
	    }
	}
	
    }

    close (F);
    O->autoflush(1);
    close(O);
}

sub formcodncod {
    my ($i1,$o1,$o2) = @_;

    open(I1,$i1);
    open(O1,">$o1");
    open(O2,">$o2");

    while(<I1>) {
	chomp;
	my @a=split;

	my $s1=substr($a[1],0,80);
	my $s2=substr($a[1],82,80);
    
	print O1 $a[0]," $s1 ",$a[2],"\n";
	print O2 $a[0]," $s2 ",$a[2],"\n";
    }
    close(I1);
    O1->autoflush(1);
    O2->autoflush(1);
    close(O1);
    close(O2);
}

sub selectfalout {

# this procedure forms a file (outfex or outfin) with a given patern 
# (e.g. t 21 -a 22, means string should have t on position 21 and 
# and not a on position 22) 
# usage: selectfalout.pl <train.falsefile> <outfalsefile> offset char1 poz1 char2 poz2 ....
# (e.g.: selectfalout.pl train.facc outfex6 44 t 21 -a 22)

    my ($trainfile,$outfile,$off,$comand)=@_;
    
    my (@nr,@consens,@poz);

    open(F,$trainfile) or die "ERROR 129: Couldn't open $trainfile for reading: $!\n";
    open(O,">$outfile") or die "ERROR 130: Couldn't open $outfile for writing: $!\n";

    while(<F>){

	my $j=0;
	my @a=split(/\s+/,$comand);

	for(my $i=0;$i<=$#a;$i+=2) {
	    $nr[$j]=0;
	    if(length($a[$i])>1) {
		$consens[$j]=substr($a[$i],1,1);
		$nr[$j]=1
		}
	    else { $consens[$j]=$a[$i];}
	
	    $poz[$j]=$a[$i+1]+$off;
	    $j++;
	}

	my $yes=1;

	my @subsir=split;

	for(my $i=0;$i<$j;$i++) {

	    # print $no[$i]," ",$consens[$i]," ",$poz[$i],"\n";

	    if($nr[$i]) {
		if(substr($subsir[1],$poz[$i],1) eq $consens[$i]) {
		    $yes=0;last;
		}
	    }
	    else {
		if(substr($subsir[1],$poz[$i],1) ne $consens[$i]) {
		    $yes=0;last;
		}
	    }
	}

	if($yes) { print O $_;}
    }
    close(F);
    O->autoflush(1);
    close(O);
}

sub formatg {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 134: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 135: Couldn't open $ex for reading: $!\n";

    my $ind=1;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=1;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start,$end}++;
	    }
	    $ind=0;
	}
    }
    
    close(E);

    open(F,$seq) or die "ERROR 136: Couldn't open $seq for reading: $!\n";

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;
	my $k;

	foreach $k (keys %isexon)
	{
	    my ($atg,$start,$end)=split(/$;/,$k);
	    if($atg eq $name)
	    {
		my $secv=substr($s,$start-13,19);
		print O "A$start $secv $name\n";
	    }
	}
    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub formfatg {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 137: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 138: Couldn't open $ex for reading: $!\n";

    my $ind=1;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_ eq "")
	{ $ind=1;}
	else   
	{
	    if($ind)
	    {
		my ($name,$start,$end)=split;
		$isexon{$name,$start}++;
	    }
	    $ind=0;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 139: Couldn't open $seq for reading: $!\n";

    my (@poz,$val);

    while(<F>)
    {
	chomp;
	my ($name,$s)=split;
	my $l=length($s);

	my $count=0;
	for(my $i=70;$i<$l-12;$i++)
	{
	    if(substr($s,$i,3) eq "atg")
	    {
		if(!$isexon{$name,$i+1}) { $poz[$count++]=$i+1;}
	    }
	}

	for(my $i=0;$i<60;$i++)
	{
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count))
	    {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }
	    
	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-13,19);
	    print O "FA$val $secv $name\n";
	    $val=0;
	}
	
    }

    close(F);
    O->autoflush(1);
    close(O);
    
}


sub formstop {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 163: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 164: Couldn't open $ex for reading: $!\n";

    my $ind=1;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_) {
	    my ($name,$start,$end)=split;
	    $isexon{$name}=$end;
	}
    }
    
    close(E);

    open(F,$seq) or die "ERROR 165: Couldn't open $seq for reading: $!\n";

    while(<F>)
    {
	chomp;

	my ($name,$s)=split;

	my $secv=substr($s,$isexon{$name}-7,19);
	print O "A",$isexon{$name}," $secv $name\n";
    }

    close(F);
    O->autoflush(1);
    close(O);

}

sub formfstop {

    my ($ex,$seq,$output)=@_;
    open(O,">$output") or die "ERROR 166: Couldn't open $output for writing: $!\n";

    open(E,$ex) or die "ERROR 167: Couldn't open $ex for reading: $!\n";

    my $ind=1;
    my %isexon;

    while(<E>)
    {
	chomp;
	if($_) {
	    my ($name,$start,$end)=split;
	    $isexon{$name}=$end;
	}
    }

    close(E);

    open(F,$seq) or die "ERROR 168: Couldn't open $seq for reading: $!\n";

    my (@poz,$val);

    while(<F>)
    {
	chomp;
	my ($name,$s)=split;
	my $l=length($s);

	my $count=0;
	for(my $i=100;$i<$l-100;$i++) {
    
	    my $stop=substr($s,$i,3);

	    if(($stop eq "taa")||($stop eq "tga")||($stop eq "tag")) {
	
		if($isexon{$name}!=$i+3) { $poz[$count++]=$i+3;}
	    }
	}

	for(my $i=0;$i<60;$i++) {
    
	    my $c=0;
	    my $j;
	    while((!$val)&&($c<2*$count)) {
		$j=int(rand($count));
		$val=$poz[$j];
		$c++;
	    }

	    $poz[$j]=0;
	
	    my $secv=substr($s,$val-7,19);
	    print O "FA$val $secv $name\n";
	    $val=0;
	}
	
    }

    close(F);
    O->autoflush(1);
    close(O);
    
}

