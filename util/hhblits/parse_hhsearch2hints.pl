#!/usr/bin/env perl

=pod

=head1 USAGE

 -infile          :s   => Output from HHblits
 -min_exons_hints :i   => Minimum number of exons with Uniprot hits before flagging sequence having a valid hit (def 2)
 -is_repeat            => HHblits was ran against a database of repeats/transposons
 -verbose              => Print out commands

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($infile, $is_repeat,$verbose);
my $min_exons_before_reporting =2;

GetOptions(
	'infile:s' => \$infile,
	'is_repeat' => \$is_repeat,
	'verbose|debug' => \$verbose,
	'min_exons_hints:i' => $min_exons_before_reporting

);


my %reporting_hash;

 my ( $homology_prob_cut, $eval_cut, $pval_cut, $score_cut,$align_col_cut, $template_aln_size_cut) = (70,1e-3,1e-6,100,50,30);
 pod2usage if !$infile || !-s $infile;
 $infile = &remove_zero_bytes( $infile) unless $infile=~/^hhr/;
 my $outfile = "$infile.results";
 die "Output already exists: $outfile\n" if -s $outfile;
 print "Post-processing $infile\n";   
 my $min_filesize = 500;
 my ( $qcounter,%hit_counter );

 die "Can't find $infile or it is too small\n" unless -s $infile && ( -s $infile ) >= $min_filesize;

 open( IN, $infile ) || die($!);

 my ($query,$id,$reverse,$start,$stop);

 while ( my $ln = <IN> ) {
  if ( $ln =~ /^\W*Query\s+(\S+)/ ) {
   $qcounter++;
   $query = $1;
   $id    = $query;
   $id =~ s/_\d+$//;
   $reverse = ($ln =~ /REVERSE SENSE/) ? 1 : 0;
   $ln =~ /\[(\d+)\s\-\s(\d+)\]/;
   $start = $1 && $1 =~ /^(\d+)$/ ? $1 : int(0);
   $stop  = $2 && $2 =~ /^(\d+)$/ ? $2 : int(0);
   next;
  }
  elsif ( $ln =~ /^\s*No Hit/ ) {
   while ( my $ln2 = <IN> ) {
    last if $ln2 =~ /^\s*$/;
    last if $ln2 =~ /^Done/;
    $ln2 =~ /^\s*(\d+)/;
    my $hit_number = $1;
    next unless $hit_number == 1;
    my ( $hit_desc, $hit_data, $hit_id );
    $hit_desc = substr( $ln2, 4, 31 );
    $hit_data = substr( $ln2, 35 );
    $hit_desc =~ /^(\S+)\s*(.*)/;
    $hit_id   = $1;
    $hit_desc = $2;

    if ($hit_desc) {
     $hit_desc =~ s/[\s\.]+$//;
     $hit_desc =~ s/\s+/_/g;
    }
    chomp($hit_data);
    $hit_data =~ s/^\s+//;
    my @data = split( /\s+/, $hit_data );
    my ( $prob, $evalue, $pvalue, $score, $structure_score, $alignment_length )
      = ( $data[0], $data[1], $data[2], $data[3], $data[4], $data[5] );
    $data[6] =~ /(\d+)\-(\d+)/;
    my $aa_start = $1;
    my $aa_stop  = $2;
    $data[7] =~ /(\d+)\-(\d+)/;
    my $hit_start = $1;
    my $hit_stop  = $2;

    if ( $data[7] =~ s/\((\d+)\)// ) {
     $data[8] = $1;
    }
    else {
     $data[8] =~ s/[\(\)]//g;
    }
    my $template_size     = $data[8];
    my $template_aln_size = abs( $hit_start - $hit_stop ) + 1;

    next if $homology_prob_cut > $prob;
    next if $eval_cut && $eval_cut < $evalue;
    next if $pval_cut && $pval_cut < $pvalue;
    next if $score_cut && $score_cut > $score;
    next if $alignment_length < $align_col_cut;
    next if $template_aln_size < $template_aln_size_cut;

    my ( $gff_start, $gff_end,$strand );
    if ( !$reverse ) {
     $gff_start = $start + ( 3 * $aa_start ) - 3;
     $gff_end   = $start + ( 3 * $aa_stop ) - 1;
     $strand = '+';
    }
    else {
     $gff_start = $start - ( 3 * $aa_stop ) + 1;
     $gff_end   = $start - ( 3 * $aa_start ) + 3;
     $strand = '-';
    }
    my $src  = $is_repeat ? 'RM'           : 'HU';
    my $type = $is_repeat ? 'nonexonpart' : 'CDSpart';
    my $prio = $is_repeat ? 6             : 5;
    my $uid  = "$hit_id.s$hit_start.e$hit_stop";
    $hit_counter{$uid}++;
    $uid .= '.n' . $hit_counter{$uid};

    my $name = $uid;
    $name .= " ($hit_desc)" if $hit_desc;
    my %hash = (
	'type' => $type,'gff_start'=>$gff_start,'gff_end'=>$gff_end,'score'=>$score,'strand'=>$strand,'src'=>$src,'prio'=>$prio,'uid'=>$uid,'hit_start'=>$hit_start,'hit_stop'=>$hit_stop,'name'=>$name,'evalue'=>$evalue
    );

    push(@{$reporting_hash{$id}{$hit_id}},\%hash);
    last;    # top hit
   }
  }
 }
 close IN;

 open( OUT,     ">$outfile" );
 open( GLIMMER, ">$outfile.glimmer" );
 open( GENEID,  ">$outfile.geneid" );
 open( GFF3,    ">$outfile.gff3" );
 open( HINTS,   ">$outfile.hints" );


 foreach my $id (sort keys %reporting_hash) {
  foreach my $hit_id (keys %{$reporting_hash{$id}}){
	my $number_of_times = scalar(@{$reporting_hash{$id}{$hit_id}});
	if ($is_repeat || $number_of_times >= $min_exons_before_reporting){
	  print OUT $id . "\n";
	  last;
	}
  }
  foreach my $hit_id (keys %{$reporting_hash{$id}}){
	my $number_of_times = scalar(@{$reporting_hash{$id}{$hit_id}});
	if ($is_repeat || $number_of_times >= $min_exons_before_reporting){
	     foreach my $hash_ref (@{$reporting_hash{$id}{$hit_id}}){
		my ($type,$gff_start,$gff_end,$score,$strand,$src,$prio,$uid,$name,$hit_start,$hit_stop,$evalue) = ( $hash_ref->{'type'},$hash_ref->{'gff_start'},$hash_ref->{'gff_end'},$hash_ref->{'score'},$hash_ref->{'strand'},$hash_ref->{'src'},$hash_ref->{'prio'},$hash_ref->{'uid'},$hash_ref->{'name'},$hash_ref->{'hit_start'},$hash_ref->{'hit_stop'},$hash_ref->{'evalue'});
		     print HINTS "$id\thhblits\t$type\t$gff_start\t$gff_end\t$score\t$strand\t.\tsrc=$src;grp=$hit_id;pri=$prio\n";
		     print GFF3  "$id\thhblits\tprotein_match\t$gff_start\t$gff_end\t$score\t$strand\t.\tID=$uid;Name=$name;Target=$hit_id $hit_start $hit_stop\n";
		     print GENEID "$id\thhblits\tsr\t$gff_start\t$gff_end\t$score\t$strand\t.\n";
	 	     if ($strand eq '-') {
			print GLIMMER "$id $gff_end $gff_start $score $evalue\n\n";
		     }else {
			print GLIMMER "$id $gff_start $gff_end $score $evalue\n\n";
		     }
	     }
	}
  }   
 }
 close OUT;
 close GLIMMER;
 close HINTS;
 close GFF3;
 close GENEID;
 if ( -s "$outfile.gff3" ) {
  system("sort -nk 4,4 $outfile.gff3| sort -s -k 1,1 > $outfile.gff3.sorted");
  rename( "$outfile.gff3.sorted", "$outfile.gff3" );
  system("sort -nk4,4 $outfile.hints|sort -s -k1,1 > $outfile.hints. ");
  rename( "$outfile.hints.", "$outfile.hints" );
 }

print "Done, see $outfile\n";

#################################################
sub remove_zero_bytes() {
 my $infile    = shift;
 my $outfile = "hhr.$infile.db";
 return $outfile if (-s $outfile);
 &process_cmd("cat $infile*.idx* > $outfile.idx");
 system("rm -f $infile*.idx*");
 &process_cmd("cat $infile* | tr -d '\\000' > $outfile");
 system("rm -f $infile*");
 return $outfile;
}

sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n" if $verbose;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

