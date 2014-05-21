#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Carp;

my $file = shift ||die;
my $delimiter = shift;


my %sort_order = (
	'First' =>1,
	'Terminal' => 1,
	'Single' =>1,
	'Internal' => 2
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

  foreach my $rec ( sort {
	# sort gene records according to scaffold and coord
      my @first  = split( "\n", $a );
      my @second = split( "\n", $b );
      my @split1 = split( "\t", $first[0] ); # gene
      my @split2 = split( "\t", $second[0] ); # gene
      return -1 if !$split2[0] || !$split2[2] || $split2[2] ne 'gene';
      return 1 if !$split1[0] || !$split1[2] || $split1[2] ne 'gene';
      return  $split1[0] cmp $split2[0] && $split1[3] <=> $split2[3];
    } @records ) {
	
	my @gene_data = split("\n",$rec);
	next if !$gene_data[2];
	@gene_data = sort {
		# sort gene lines according to sort_order
	      my @split1 = split( "\t", $a );
	      my @split2 = split( "\t", $b );
	      return -1 if !$split2[0] || !$split2[2] || !$split2[3] || !$sort_order{$split2[2]};
	      return 1 if !$split1[0] || !$split1[2] || !$split1[3] || !$sort_order{$split1[2]};
	      return $sort_order{$split1[2]} <=> $sort_order{$split2[2]} && $split1[3] <=> $split2[3];
	}(@gene_data);
	$rec = join("\n",@gene_data);
    print OUT $rec . "\n".$/;
  }
  close GFF;
  close OUT;
  $/ = $orig_sep;
  system("uniq $gff.sorted > $gff.fixed");
  unlink("$gff.sorted");
#  rename( "$gff.sorted", $gff );
}

sub get_gff_delimiter() {
 my $delimiter;
 my $file = shift;
 return if !-s $file;
 open( IN, $file ) || confess("Can't find file $file");
 my $skip = <IN>;
 while ( my $ln = <IN> ) {
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
