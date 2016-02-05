#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
my @files = @ARGV;

die unless $files[1] && -s $files[1];

#  $cfg{$source}->{ $feature_headers[$i] }->{'value'}  = $score

my (%common,%sources,%srcparam,$manual_str,%cfg,$common_print);
my @feature_order =
  qw/start stop tss tts ass dss exonpart exon intronpart intron CDSpart CDS UTRpart UTR irpart nonexonpart genicpart/;

foreach my $f (@files){
	&parse_cfg($f);
}

$common_print='[SOURCES]'."\nM";
foreach my $src (keys %sources){
	$common_print.= ' '.$src;
}


$common_print.= "\n\n";
$common_print.= '[SOURCE-PARAMETERS]'."\n";
foreach my $ln (keys %srcparam){
	$common_print.= $ln;
}


$common_print.= "\n";
foreach my $ln (keys %common){
	$common_print.= $ln;
}

my $number_of_cfgs_to_print = 1;
my (%prints);
foreach my $feature (@feature_order){
	foreach my $bonus (1,keys %{$cfg{$feature}{'bonus'}}){
		foreach my $malus (1,keys %{$cfg{$feature}{'malus'}}){
			foreach my $lc_malus (1,keys %{$cfg{$feature}{'local_malus'}}){
				my $p = "$feature $bonus $malus $lc_malus ".$manual_str;
				foreach my $src (keys %sources){
					if (exists($cfg{$feature}{$src}) && scalar(keys %{$cfg{$feature}{$src}})>1){
						foreach my $value (1,keys %{$cfg{$feature}{$src}}){
							$p .= " $src 1 $value";
						}
					}else{
						my @values = keys %{$cfg{$feature}{$src}};
						$values[0] = 1 unless $values[0];
						$p .= " $src 1 $values[0]";
					}
				}
				push(@{$prints{$feature}},$p ."\n");
			}
		}
	}
	$number_of_cfgs_to_print = $number_of_cfgs_to_print * scalar(@{$prints{$feature}});
}

my (@fhs);
my %cfg_vals;
my %feat_tupple;
for (my $cfg_i=int(0);$cfg_i<$number_of_cfgs_to_print ; $cfg_i++){
	my $outfile = "out.cfg.".$cfg_i;
	open (my $fh,">$outfile");
	$fhs[$cfg_i] = $fh;
	print $fh $common_print;
	foreach my $feature (@feature_order){
		my @values = @{$prints{$feature}};
		if (scalar(@values)==1){
			print $fh $values[0];
		}else{
			@{$cfg_vals{$feature}} = @values;
			$feat_tupple{$feature}=1;
		}
	}	
}

my $num_feat_tupple =  scalar(keys %feat_tupple);
warn $num_feat_tupple;
warn Dumper \%cfg_vals;

my %print2;
for (my $cfg_i=int(0);$cfg_i<$number_of_cfgs_to_print ; $cfg_i++){
	my $fh = $fhs[$cfg_i];
	for (my $f=0;$f<$num_feat_tupple;$f++){
		my $feature = ${keys %feat_tupple}[$f];
die $feature;
		my @Fvalues = @{$cfg_vals{$feat_tupple{$feature}}};
		print $fh $Fvalues[$cfg_i] if $Fvalues[$cfg_i];
	}
}
for (my $cfg_i=$num_feat_tupple;$cfg_i<$number_of_cfgs_to_print ; $cfg_i++){
	my $fh = $fhs[$cfg_i];
	for (my $f=0;$f<$num_feat_tupple;$f++){
		my @Fvalues = @{$cfg_vals{$feat_tupple{$f}}};
		print $fh $Fvalues[$cfg_i] if $Fvalues[$cfg_i];
	}
}


foreach my $fh (@fhs){
	close $fh;
}

#die Dumper \%prints;




######
sub parse_cfg(){
my ($file) = shift;
open (my $fh,$file);
while (my $ln=<$fh>){
	next if $ln=~/^\s*$/;
	if ($ln=~/^\[SOURCES\]/){
		$ln = <$fh>;chomp($ln);
		my @sources_k = split(/\s+/,$ln);
		foreach my $s (@sources_k){
			$sources{$s}++ unless $s eq 'M';
		}
	}

	if ($ln=~/^#/){
		$common{$ln}++;
		next;
	}
	if ($ln=~/\[SOURCE-PARAMETERS\]/){
		while ($ln!~/^\s*$/){
			$ln=<$fh>;
			$srcparam{$ln}++;
		}
	}
	if ($ln=~/^\[GENERAL\]/){
		while ($ln!~/^\s*$/){
			$ln=<$fh>;
			chomp($ln);
			my ($feature,$bonus,$malus,$local_malus,@data) = split(/\s+/,$ln);
			#0-2 are always the same
			#3-5 are the source's
			next unless scalar(@data) > 5;
			my @manual_data = splice(@data,0,3);
			$manual_str = join(' ',@manual_data) unless $manual_str;
			$cfg{$feature}{'bonus'}{$bonus}++ unless $bonus == 1;
			$cfg{$feature}{'malus'}{$malus}++ unless $malus == 1;
			$cfg{$feature}{'local_malus'}{$local_malus}++ unless $local_malus == 1;
			$cfg{$feature}{$data[0]}{$data[2]}++ unless $data[2] == 1;
		}

	}
}

close $fh;
}

