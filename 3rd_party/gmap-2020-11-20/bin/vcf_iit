#! /usr/bin/perl

use IO::File;
use Getopt::Std;
undef $opt_C;			# Include only those lines that are marked COMMON
getopts("Cv:n:");

if (defined($opt_n)) {
    # Names file in 2-column format, which should be the same as for gmap_build
    $FP = new IO::File($opt_n) or die "Cannot open $opt_n";
    while (defined($line = <$FP>)) {
	chop $line;
	($old,$new) = split /\s+/,$line;
	if (defined($old) && defined($new) && $new =~ /\S/) {
	    $contig_newname{$old} = $new;
	}
    }
    close($FP);
}


while (defined($line = <>)) {
      if ($line =~ /^\#/) {
	  # Skip
      } else {
	  chop $line;
	  @fields = split /\t/,$line;

	  $chr = $fields[0];
	  if (defined($contig_newname{$chr})) {
	      $chr = $contig_newname{$chr};
	  }

	  $chrpos = $fields[1];
	  $id = $fields[2];
	  $ref_allele = $fields[3];
	  $alt_allele = $fields[4];

	  $version = get_info("dbSNPBuildID",$fields[7]);
	  if (defined($opt_v) && $version > $opt_v) {
	      $wantp = 0;
#	  } elsif (defined(get_info("VLD",$fields[7]))) {
#	      # Validated
#	      $wantp = 1;
#	  } elsif (defined(get_info("SLO",$fields[7]))) {
#	      # Submitter link-out
#	      $wantp = 1;
#	  } elsif (defined(get_info("GNO",$fields[7]))) {
#	      # Individual genotype
#	      $wantp = 0;
#	  } else {
#	      $wantp = 1;

	  } elsif (!defined($opt_C)) {
	      $wantp = 1;
	  } elsif (defined(get_info("COMMON",$fields[7]))) {
	      $wantp = 1;
	  } else {
	      $wantp = 0;
	  }

	  if ($wantp == 1 && $ref_allele =~ /^[ACGT]$/ && $alt_allele =~ /^[ACGT]$/) {
	      if ($ref_allele lt $alt_allele) {
		  print ">$id $chr:$chrpos ";
		  print $ref_allele . $alt_allele;
		  print "\n";
	      } elsif ($alt_allele lt $ref_allele) {
		  print ">$id $chr:$chrpos ";
		  print $alt_allele . $ref_allele;
		  print "\n";
	      } else {
		  print STDERR "Ref allele $ref_allele and alt allele $alt_allele for $id at $chr:$chrpos are the same\n";
	      }
	  }
      }
}

exit;


sub get_info {
    my ($desired_key, $info) = @_;

    foreach $binding (split ";",$info) {
	if ($binding =~ /(\S+)=(\S+)/) {
	    $key = $1;
	    $value = $2;
	    if ($key eq $desired_key) {
		return $value;
	    }
	} else {
	    $key = $binding;
	    if ($key eq $desired_key) {
		return 1;
	    }
	}
    }
    return;
}

exit;

