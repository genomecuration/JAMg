#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Carp;

my $file = shift ||die;
my $delimiter = shift;


my %sort_order = (
	'gene' =>1,
	'mRNA' =>2,
	'three_prime_UTR' => 3,
	'five_prime_UTR' => 3,
	'exon' =>3,
	'intron' =>3,
	'CDS' =>5,
	'tss' => 3,
	'tts' => 3,
);

$delimiter =&get_gff_delimiter($file) unless $delimiter;
print "Using delimiter $delimiter\n";

&sort_gff3($file,$delimiter);

sub sort_gff3() {
  my $gff       = shift;
  my $delimiter = shift;
  my $orig_sep  = $/;
  open( GFF, $gff )           || die;
  open( OUT, ">$gff.sorted" ) || die;
  $/ = $delimiter;
  my @records = <GFF>;
  chomp(@records);
  print "Processing ".scalar(@records)." records\n";
  foreach my $rec ( sort {
	# sort gene records according to scaffold and coord
      my @first  = split( "\n", $a );
      my @second = split( "\n", $b );
      return 1 if !$first[0];
      return -1 if !$second[0];
      my @split1 = split( "\t", $first[0] ); # gene
      my @split2 = split( "\t", $second[0] ); # gene
      return -1 if !$split2[0] || !$split2[2] || $split2[2] ne 'gene';
      return 1 if !$split1[0] || !$split1[2] || $split1[2] ne 'gene';
      return  $split1[0] cmp $split2[0] || $split1[3] <=> $split2[3];
    } @records ) {
	
	my @gene_data = split("\n",$rec);
	my @n;
	for (my $i=0;$i<@gene_data;$i++){
		my @d = split("\t",$gene_data[$i]);
		if ($d[2] && $d[2] eq 'transcript'){
			$d[2] = 'mRNA';
			$gene_data[$i] = join("\t",@d);
		}
		if ($d[2] && $d[2] eq 'transcription_start_site'){
			$d[2] = 'tss';
			$gene_data[$i] = join("\t",@d);
		}
		if ($d[2] && $d[2] eq 'transcription_end_site'){
			$d[2] = 'tts';
			$gene_data[$i] = join("\t",@d);
		}
		push(@n,$gene_data[$i]) if $d[2] && $sort_order{$d[2]};
	}
	@gene_data = @n;
	next if !$gene_data[2];	# at least one gene, mrna and exon
	@gene_data = sort {
		# sort gene lines according to sort_order
	      my @split1 = split( "\t", $a );
	      my @split2 = split( "\t", $b );
	      return -1 if !$split2[0] || !$split2[2] || !$split2[3] || !$sort_order{$split2[2]};
	      return 1 if !$split1[0] || !$split1[2] || !$split1[3] || !$sort_order{$split1[2]};
	      return $sort_order{$split1[2]} <=> $sort_order{$split2[2]} || $split1[3] <=> $split2[3];
	}(@gene_data);
	$rec = join("\n",@gene_data);
    print OUT $rec . "\n\n";
  }
  close GFF;
  close OUT;
  $/ = $orig_sep;
  system("uniq $gff.sorted > $gff.fixed");
  rename( "$gff.fixed", "$gff.sorted" );
}

sub get_gff_delimiter() {
 my $delimiter;
 my $file = shift;
 return if !-s $file;
 open( IN, $file ) || confess("Can't find file $file");
 my $skip = <IN>;
 while ( my $ln = <IN> ) {
  last if $ln =~ /^##FASTA/;
  if ( $ln =~ /^\s*$/ ) {

   # note that delimiter of an empty line is actually two \n\n
   $delimiter = $ln . "\n";
   last;
  }
  elsif ( $ln =~ /^#/ ) {
   $delimiter = $ln;
   last;
  }
 }
 close IN;
 confess "I don't know what delimiter to use for $file" if !$delimiter;
 return $delimiter;
}
