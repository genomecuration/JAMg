#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
my @files = @ARGV;

die unless $files[1] && -s $files[1];

#  $cfg{$source}->{ $feature_headers[$i] }->{'value'}  = $score

my (%common,%sources,%srcparam,$manual_str,%cfg);
my @feature_order =
  qw/start stop tss tts ass dss exonpart exon intronpart intron CDSpart CDS UTRpart UTR irpart nonexonpart genicpart/;

foreach my $f (@files){
	&parse_cfg($f);
}

open (OUT,">out.cfg");
print OUT '[SOURCES]'."\nM";
foreach my $src (keys %sources){
	print OUT ' '.$src;
}


print OUT "\n\n";
print OUT '[SOURCE-PARAMETERS]'."\n";
foreach my $ln (keys %srcparam){
	print OUT $ln;
}


print OUT "\n";
foreach my $ln (keys %common){
	print OUT $ln;
}

my (%prints);
foreach my $feature (@feature_order){
	print OUT $feature;
	foreach my $bonus (1,keys %{$cfg{$feature}{'bonus'}}){
		foreach my $malus (1,keys %{$cfg{$feature}{'malus'}}){
			foreach my $lc_malus (1,keys %{$cfg{$feature}{'local_malus'}}){
				my $p;
				$p = "$feature $bonus $malus $lc_malus ".$manual_str;
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
}


die Dumper \%prints;

close OUT;



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

