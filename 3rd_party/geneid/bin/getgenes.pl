#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %switches;
my $nonred=0;
my $onlynonred=0;

getopts("hrn",\%switches);

if (exists($switches{h})) {
  print "\ngetgenes.pl [-hr] <path2goldenpath> <genecoordinates>\n\n";
  print "synopsis: retrieves genes in forward from a goldenpath genome\n";
  print "          by reading gene coordinates in goldenpath format\n\n";
  print "parameters: <path2goldenpath> path where the goldenpath genome is\n";
  print "            <genecoordinates> file with genes coordinates in goldenpath format\n";
  print "\noptions: -h this help\n";
  print "         -r annotate redundant genes, i.e. genes with non-unique\n";
  print "            names and/or transcript coordinates that overlap with\n";
  print "            at least one other gene\n";
  print "         -n do not print redundant genes (only works with -r)\n\n";

  exit(0);
}

if (exists($switches{r})) {
  $nonred=1;
}

if (exists($switches{n})) {
  if ($nonred) {
    $onlynonred=1;
  } else {
    print "getgenes.pl (error): -n only works with -r\n";
    exit(1);
  }
}

if (scalar(@ARGV) < 2) {
  print "getgenes.pl [-hrn] <path2goldenpath> <genecoordinates>\n";
  exit(1);
}

my $path2gpath = $ARGV[0];
my $genesfname = $ARGV[1];
my $prevgene = "x";
my $prevchro = "x";
my $trail = "";
my %txcoord;
my %genenames;

chomp($path2gpath);
chomp($genesfname);

if (!open(REFGENE,"< $genesfname")) {
  print "getgenes.pl: impossible to open $genesfname\n";
  exit(1);
}

if ($nonred) {
  if (!open(REDGENE,"> redgenes.txt")) {
    print "getgenes.pl: impossible to create redgenes.txt\n";
    exit(1);
  }
}

while (<REFGENE>) {

  m/([\w\-\.:]+)\s+([\w\.\-:]+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([^\n]+)/;

  my $name = $1;
  my $chro = $2;
  my $stra = $3;
  my $txSt = $4;
  my $txEn = $5;
  my $cdsS = $6;
  my $cdsE = $7;
  my $exoC = $8;
  my @exon = ($9 =~ m/(\d+)/g);
  my $cdsLe = $cdsE - $cdsS;
  my $cdsoffset = $cdsS - $txSt;
  my $redundant=0;
  my $i = 0; # exon counter
  my $j = 0; # intron counter

  if ($nonred) {
    my $origname;
    my $origchro;
    my @origexon;
    my $overlap=0;

    if (geneoverlaps($name,$chro,$txSt,$txEn,\@exon,\$origname,\@origexon)) {
      print REDGENE "$name\t$origname\t$chro";
      my $iorig=0;
      my $origexoC=scalar(@origexon)/2;
      $i=0;

      while ($iorig < $origexoC && $i < $exoC) {
        if ($exon[$i]!=$origexon[$iorig] && $exon[$i+$exoC]!=$origexon[$iorig+$origexoC]) {
          print REDGENE "\t$exon[$i],$exon[$i+$exoC]";
          if ($exon[$i] < $origexon[$iorig]) {
            $i = $i + 1;
          } else {
            $iorig = $iorig + 1;
          }
        } elsif ($exon[$i]!=$origexon[$iorig]) {
          print REDGENE "\t$exon[$i],-";
          $iorig = $iorig + 1;
          $i = $i + 1;
        } elsif ($exon[$i+$exoC]!=$origexon[$iorig+$origexoC]) {
          print REDGENE "\t-,$exon[$i+$exoC]";
          $iorig = $iorig + 1;
          $i = $i + 1;
        } else {
          $iorig = $iorig + 1;
          $i = $i + 1;
        }
      }

      while ($i < $exoC) {
        print REDGENE "\t$exon[$i],$exon[$i+$exoC]";
        $i = $i + 1;
      }
      print REDGENE "\n";
      $overlap=1;
    }

    if (exists($genenames{$name})) {
      my $searchname=$name;
      $genenames{$name} = $genenames{$name} + 1;
      if ($overlap) {
        $name = $name . ".$origname"; 
      }
      $name = $name . ".$genenames{$searchname}";
      $redundant = 1;
    } else {
      $genenames{$name} = 0;
      if ($overlap) {
        $name = $name . ".$origname"; 
        $redundant = 1;
      }
    }
  }

  if (!$onlynonred || ($onlynonred && !$redundant)) {

    my $call = "./chromosomechunk ".$path2gpath.$chro." ".$txSt." ".($txEn-$txSt);
if ($txSt==0) { print STDERR "SOMETHINGWRONG with $name: $call\n"; exit(1);}
    my $genomic = `$call`;
    my $genomicLe = length($genomic);
    my $cdseq = "";

    if ($genomicLe == 0) {
      print STDERR "getgenes.pl: gene of 0 length ($name), $call\n";
      next;
    }

    if ($genomicLe != $txEn - $txSt) {
      print STDERR "getgenes.pl: length mismatch ($name)\n";
      next;
    }

    for ($i=0;$i<$exoC;$i++) {

      my $utrB = 0;
      my $utrA = 0;
      my $utrS = 0;
      my $utrL = 0;
      my $exSt = $exon[$i] - $cdsS;
      my $exLe = $exon[$i+$exoC] - $exon[$i];
      my $exTy = "Internal";

      if ($exSt+$exLe > 0 && $exSt < $cdsLe) { # cds
  
        if ($exSt <= 0 || $i == 0) {
          if ($stra eq '+') {
            $exTy = "First";
          } else {
            $exTy = "Terminal";
          }
        }

        if ($exSt+$exLe >= $cdsLe || $i == $exoC-1) {
          if ($stra eq '+') {
            $exTy = "Terminal";
          } else {
            $exTy = "First";
          }
        }

        if ($exSt <= 0 && $exSt+$exLe >= $cdsLe) {
          $exTy = "Single";
        }

        if ($exSt < 0) {
          $utrB = 1;
          $utrS = $exSt;
          $utrL = abs($exSt);
          $exLe = $exLe - abs($exSt);
          $exSt = 0;
        }

        if ($exSt+$exLe > $cdsLe) {
          $utrA = 1;
          $utrS = $cdsLe;
          $utrL = $exLe - ($cdsLe - $exSt);
          $exLe = $cdsLe - $exSt;
        }

        my $iex;
        my $seq = substr($genomic,$exSt+$cdsoffset,$exLe);

        $seq = lc($seq);

        if ($stra eq '+') {  # forward

          if ($utrB) {
            my $iutr = $i+1;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
          }
  
          $iex = $i+1;
          $cdseq = $cdseq . "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i]\t$exon[$i+$exoC]\n";

          if ($utrA) {
            my $iutr = $i+1;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
          }

        } else {             # reverse

          if ($utrB) {
            my $iutr = $exoC-$i;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $utrs =~ tr/acgt/tgca/;
            $utrs = reverse($utrs);
            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n" . $cdseq;
          }

          $iex = $exoC-$i;
          $seq =~ tr/acgt/tgca/;
          $seq = reverse($seq);
          $cdseq = "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i+$exoC]\t$exon[$i]\n" . $cdseq;

          if ($utrA) {
            my $iutr = $exoC-$i;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);

            $utrs = lc($utrs);
            $utrs =~ tr/acgt/tgca/;
            $utrs = reverse($utrs);
            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n" . $cdseq;
          }

        }

        if ($exTy ne "Single" && (($exTy ne "Terminal" && $stra eq '+') ||
           ($exTy ne "First" && $stra eq '-'))) {

          my $inSt = $exon[$i+$exoC] - $cdsS;
          my $inLe = $exon[$i+1] - $exon[$i+$exoC];

          if ($inSt+$inLe > 0 && $inSt < $cdsLe) {

            if ($inSt < 0) {
              print "getgenes.pl: intron out of range! (1)\n";
              exit(1);
            }

            if ($inSt+$inLe > $cdsLe) {
              print "getgenes.pl: intron out of range! (2)\n";
              exit(1);
            }

            $seq = substr($genomic,$inSt+$cdsoffset,$inLe);
            $seq = "\L$seq";

            my $iIn;

            if ($stra eq '+') { # forward
              $iIn = $j+1;
              $cdseq = $cdseq . "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n";
            } else {
              $iIn = $exoC-$j-1;
              $seq =~ tr/acgt/tgca/;
              $seq = reverse($seq);
              $cdseq = "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n" . $cdseq;
            }

          } else {
            print "getgenes.pl: intron out of range! (3)\n";
            if ($inSt+$inLe <= 0) {
              print "getgenes.pl: intron in 5' UTR\n";
            } else {
              print "getgenes.pl: intron in 3' UTR\n";
            }

            exit(1);
          }

          $j = $j + 1;
        }

      } else {  # UTRs
  
        $exSt = $exon[$i] - $txSt;
        $exLe = $exon[$i+$exoC] - $exon[$i];
  
        my $utrs = substr($genomic,$exSt,$exLe);
  
        if ($stra eq '+') {  # forward
          my $iutr = $i+1;

          $utrs = lc($utrs);
          $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n";
        } else {             # reverse
          my $iutr = $exoC-$i;

          $utrs = lc($utrs);
          $utrs =~ tr/acgt/tgca/;
          $utrs = reverse($utrs);
          $cdseq = "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n" . $cdseq;
        }

      }
    }

    print $cdseq;

  } elsif ($onlynonred) {
    print STDERR "$name\n";
  }

}

if ($nonred) {
  close(REDGENE);
}

# we use a binary tree to store the gene transcript coordinates and check whether a given
# gene overlaps with the rest efficiently

sub geneoverlaps {
  my $name = $_[0];
  my $chro = $_[1];
  my $txSt = $_[2];
  my $txEn = $_[3];
  my @exon = @{$_[4]};

  if (!exists($txcoord{$chro})) {
    $txcoord{$chro} = [ undef, undef, $txSt, $txEn, $name, [@exon] ];
    return 0;
  }

  my $p=$txcoord{$chro};
  my $prev=undef;

  while ($p) {
    if (($txSt >= $p->[2] && $txSt <= $p->[3]) || ($txEn >= $p->[2] && $txEn <= $p->[3])) {
      ${$_[5]} = $p->[4];
      @{$_[6]} = @{$p->[5]};
      return 1;
    }

    $prev = $p;
    if ($txSt <= $p->[2]) {
      $p = $p->[0];
    } else {
      $p = $p->[1];
    }
  }

  if (!$prev) {
    print STDERR "getgenes.pl: something is wrong in the binary tree\n";
    exit(1);
  }

  if ($txSt <= $prev->[2]) {                             # insert
    $prev->[0] = [ undef, undef, $txSt, $txEn, $name, [@exon] ];
  } else {
    $prev->[1] = [ undef, undef, $txSt, $txEn, $name, [@exon] ];
  }

  return 0;
}
