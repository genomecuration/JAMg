#!/usr/bin/perl -w

# This script validates splice sites across introns as annotated in GTF file
# That is, it tracks down splice site pairs like GT-AC or AT-AG that contain
# valid splice sites which are paired inconsistently. This script is intended
# to be used after validate_gtf.pl in creating training and evaluation sets
# to augment bad genes list produced by validate_gtf.pl

use strict;
use lib '.';
use GTF;

my $usage = "usage: $0 <gtf file> <sequence file> <bad genes list>\n";
die $usage unless (@ARGV == 3);
my ($filename,$seqfile,$badlist) = @ARGV;

my %compl = ("A" => "T", "T" => "A", "G" => "C", "C" => "G");

my $seq = "";
open(IN, "$seqfile");
while (<IN>)
  {
     chop $_;
     next if (/^>/);
     $seq .= uc($_);
  }
close(IN);

my $gtf = GTF::new({gtf_filename => $filename,
		    warning_fh   => \*STDERR});
my $genes = $gtf->genes;
my $gene;

my %bad = ();
my @bad_genes = ();
open(IN, "$badlist");
while (<IN>)
  {
    chop $_;
    if (/^(NM_\S+\.a)$/)
      {
        push(@bad_genes, $1);
        $bad{$1} = 1;
      }
  }
close(IN);

foreach $gene (@$genes)
  {
    my $transcripts = $gene->transcripts;
    my $tx;

    foreach $tx (@$transcripts)
      {
        next if (defined $bad{$tx->id});

        my $ps = ($tx->strand eq "+");
        my $exons = $tx->cds;
        my $utr3 = $tx->utr3;
        my $utr5 = $tx->utr5;
        my $stops = $tx->stop_codons;
        my $starts = $tx->start_codons;
        my $p_exons;
        
        foreach $p_exons ($utr5, $starts, $exons, $stops, $utr3)
          {
            my $i;
            for ($i = 1; $i <= $#$p_exons; $i++)
               {
                 next unless ($$p_exons[$i]->start - $$p_exons[$i-1]->stop > 4);

                 my $splice5p = ($ps) ? substr($seq, ($$p_exons[$i-1]->stop), 2)
                                      : substr($seq, ($$p_exons[$i]->start-3), 2);

                 my $splice3p = ($ps) ? substr($seq, ($$p_exons[$i]->start-3), 2)
                                      : substr($seq, ($$p_exons[$i-1]->stop), 2);

                 &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                    ($gene->gene_id), ($tx->id));
               }
          }

        # this assumes the innermost start codon is included in CDS
        # therefore, no checks are run on CDS - start_codon junction
        # also, the check is done only on introns that are at least 
        # 4 bp long (at least 2 bp per each splice site)
        if ($ps)
          {
            if (($#$utr5 > -1) && ($#$starts > -1))
              {
                 if (($$starts[0]->start - $$utr5[$#$utr5]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$utr5[$#$utr5]->stop), 2);
                      my $splice3p = substr($seq, ($$starts[0]->start-3), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
            if (($#$exons > -1) && ($#$stops > -1))
              {
                 if (($$stops[0]->start - $$exons[$#$exons]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$exons[$#$exons]->stop), 2);
                      my $splice3p = substr($seq, ($$stops[0]->start-3), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
            if (($#$stops > -1) && ($#$utr3 > -1))
              {
                 if (($$utr3[0]->start - $$stops[$#$stops]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$stops[$#$stops]->stop), 2);
                      my $splice3p = substr($seq, ($$utr3[0]->start-3), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
          }
        else
          {
            if (($#$utr3 > -1) && ($#$stops > -1))
              {
                 if (($$stops[0]->start - $$utr3[$#$utr3]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$stops[0]->start-3), 2);
                      my $splice3p = substr($seq, ($$utr3[$#$utr3]->stop), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
            if (($#$stops > -1) && ($#$exons > -1))
              {
                 if (($$exons[0]->start - $$stops[$#$stops]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$exons[0]->start-3), 2);
                      my $splice3p = substr($seq, ($$stops[$#$stops]->stop), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
            if (($#$starts > -1) && ($#$utr5 > -1))
              {
                 if (($$utr5[0]->start - $$starts[$#$starts]->stop) > 4)
                   {
                      my $splice5p = substr($seq, ($$utr5[0]->start-3), 2);
                      my $splice3p = substr($seq, ($$starts[$#$starts]->stop), 2);
                      &check_splice_pair($splice5p, $splice3p, $ps, \@bad_genes, 
                                         ($gene->gene_id), ($tx->id));
                   }
              }
          }      
      }
  }

my $bad_gene;
foreach $bad_gene (@bad_genes)
  {
    print "$bad_gene\n";
  }

exit;

sub check_splice_pair
{
  my($s5p, $s3p, $is_pos, $bad, $gene_id, $tx_id) = @_;

  if (!$is_pos)
    {
      $s5p = &rev_comp($s5p);
      $s3p = &rev_comp($s3p);
    }

  if (($s5p eq "GT") || ($s5p eq "GC"))
    {
       if ($s3p ne "AG")
         {
           print "Inconsistent splice site pair detected: $s5p - $s3p\n";
           print "Strand: ".(($is_pos) ? "+" : "-")." Gene :".$gene_id.
                 ", transcript: ".$tx_id."\n";
           push(@$bad, $tx_id);
         }
    }
  elsif ($s5p eq "AT")
    {
       if ($s3p ne "AC")
         {
           print "Inconsistent splice site pair detected: $s5p - $s3p\n";
           print "Strand: ".(($is_pos) ? "+" : "-")." Gene: ".$gene_id.
                 ", transcript: ".$tx_id."\n";
           push(@$bad, $tx_id); 
         }
    }
  else
    {
       print "Invalid splice site pair detected: $s5p - $s3p\n";
       print "Strand: ".(($is_pos) ? "+" : "-")." Gene :".$gene_id.
             ", transcript: ".$tx_id."\n";
       push(@$bad, $tx_id);
    }
}

sub rev_comp
{
  my($splice) = uc($_[0]);
  my($rc_splice, $i);

  for ($i = 1; $i > -1; $i--)
    {
      my($nucl) = substr($splice, $i, 1);
      if ($nucl =~ /^[ATGC]$/)
        {
          $rc_splice .= $compl{$nucl};
        }
      else
        {
          $rc_splice .= "N";
        }
    }

  return $rc_splice;
}
__END__
