#!/usr/bin/env perl

=pod

=USAGE

        'fill_in=s' => 
        'longreads=s{1,}' => 
        'min_gap:i' => 
        'max_gap:i' => 
        'capture_length:i' => 
 

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;

my ($fill_in,$longreads,$do_help) ;
my $gap_size_min = 1;
my $gap_size_max = 1e4;
my $capture_length = 1000;
my $identity_cutoff = 98;
my $cpus = 1;
GetOptions ( 
	'h|?' => \$do_help,
	'fill_in=s' => \$fill_in,
	'longreads=s' => \$longreads,
	'min_gap:i' => \$gap_size_min,
	'max_gap:i' => \$gap_size_max,
	'capture_length:i' => \$capture_length,
        'identity_min:f' => \$identity_cutoff,
	'cpus:i'	=>\$cpus
) || pod2usage(2);
pod2usage(1) unless !$do_help || ($fill_in && -s $fill_in && $longreads && -s $longreads);

my ($lastz_db_exec,$lastz_aln_exec ) = &check_program('lastdb','lastal');

#split fill_in
my $lastdb = $longreads.'.contigdb';
&process_cmd($lastz_db_exec ." -c $lastdb $longreads") unless -s "$lastdb.bck";
my $query_read_file = &split_gap_fasta($fill_in);

my $lastz_aln_opts = ' -S 1 -T 1 -f BlastTab -l 200 -P '.$cpus;
my $lastz_output = $query_read_file.'_vs_'.basename($longreads).'.lastz';
&process_cmd($lastz_aln_exec ." $lastz_aln_opts $lastdb $query_read_file > $lastz_output") unless -s $lastz_output;

&process_lastz_blast($lastz_output);


#########################
sub parse_fasta_to_hash(){
  my $fasta = shift;
  die unless $fasta && -s $fasta;
  my $orig_sep = $/;
  $/ = '>';
  my %hash;
  open( IN, $fasta ) || die( "Cannot open $fasta : " . $! );
  while ( my $record = <IN> ) {
	chomp($record);
	next unless $record;
	my @lines = split( $orig_sep, $record );
	my $id    = shift(@lines);
	$id  =~s/\s.+//g;
	my $seq   = join( '', @lines );
	$seq =~s/\s+$//;
	next unless $id && $seq;
	$hash{$id} = $seq;
  }
 $/ = $orig_sep;
 return \%hash;
}

sub process_lastz_blast(){
  my ($file) = @_;
  my $output = $file.".parsed";
  my (%hash1,%hash2);

  open (IN,$file);
  while (my $ln=<IN>){
	next unless $ln && $ln!~/^#/;
	my @data = split("\t",$ln);
	# query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	#      0       1            2            3                    4      5         6          7         8      9      10      11
	next unless $data[11];
	chomp($data[11]);	
	next unless $data[2] >= $identity_cutoff || $data[3] >= $capture_length * 0.8;
	# scaffold:1:A:2491:3491:1270
        #0  'scaffold', # scaffold id
        #1  '1',    # subseq count
        #2  'A',    # pair, A or B
        #3  '2491', # start 1-coord
        #4  '3491', # end 1-coord
        #5  '1270'  # gap size between A & B
	# get rid of bloody e notation
        $data[11] = sprintf("%.10g", $data[11]);
	my @parsed_query = split(':',$data[0]);

	# need to improve the hit parsing here to ensure we get something meaningful in the face of repeats
	# go through and see if there is a pair hitting same subject
	# if there is more than one, pick the one with highest alignment length, identity, score
	

	if (    
		!$hash1{$data[1]}{$parsed_query[0].':'.$parsed_query[1]}{$parsed_query[2]} ||
		( $data[2] > $hash1{$data[1]}{$parsed_query[0].':'.$parsed_query[1]}{$parsed_query[2]}{'identity'} 
		&& $data[11] > $hash1{$data[1]}{$parsed_query[0].':'.$parsed_query[1]}{$parsed_query[2]}{'score'} )
		){
		# query can be reversed
		my ($q_start,$q_end,$orient) = $data[6] < $data[7] ? ($data[6],$data[7],'1') : ($data[7],$data[6],'-1');

		$hash1{$data[1]}{$parsed_query[0]}{$parsed_query[1]}{$parsed_query[2]} = 
		 {
			'alignment_length' => int($data[3]),
			'scaffold' => $parsed_query[0],
			'subseq' => int($parsed_query[1]),
			'fragside' => $parsed_query[2],
			'identity' => $data[2],
			'score' => int($data[11]),
			'q_real_start_1' => $q_start + $parsed_query[3] - 1, 
			'q_real_end_1' => $q_end + $parsed_query[3] , # length aligned hit, remove the rest?
			'subject' => $data[1],
			's_real_start_1' => int($data[8]), 
			's_real_end_1' => int($data[9]),
			'orientation' => $orient,
			'orig_gap_length' => int($parsed_query[5])
		 };
		# this is the subject seq, fyi
		# my $subject_seq = substr($longread_hash->{$data[1]}, $data[8] -1 , ($data[9] - $data[8] + 1 )  );
		#die $subject_seq;
	}
  }
  close IN;
#die Dumper \%hash1;
 # each subject must get both subseq fragments. but each fragment pair must be linked to the same and only one subject
  foreach my $subject (keys %hash1){
	foreach my $scaffold (keys %{$hash1{$subject}}){
	   foreach my $subseq (keys %{$hash1{$subject}{$scaffold}}){
		next unless $hash1{$subject}{$scaffold}{$subseq}{'A'} && $hash1{$subject}{$scaffold}{$subseq}{'B'} 
                     && $hash1{$subject}{$scaffold}{$subseq}{'A'}{'orientation'} == $hash1{$subject}{$scaffold}{$subseq}{'B'}{'orientation'};
		# if there is more than one, pick the one with highest alignment length, identity, score
		my ($new_aln_length,$new_score) = (
			$hash1{$subject}{$scaffold}{$subseq}{'A'}{'alignment_length'} + $hash1{$subject}{$scaffold}{$subseq}{'B'}{'alignment_length'},
			$hash1{$subject}{$scaffold}{$subseq}{'A'}{'score'} + $hash1{$subject}{$scaffold}{$subseq}{'B'}{'score'},
		   );
		if (!$hash2{$scaffold}{$subseq}{'A'}){
			$hash2{$scaffold}{$subseq}{'A'} = $hash1{$subject}{$scaffold}{$subseq}{'A'};
			$hash2{$scaffold}{$subseq}{'B'} = $hash1{$subject}{$scaffold}{$subseq}{'B'};
		}else{
		  my ($current_aln_length,$current_score) = (
			$hash2{$scaffold}{$subseq}{'A'}{'alignment_length'} + $hash2{$scaffold}{$subseq}{'B'}{'alignment_length'},
			$hash2{$scaffold}{$subseq}{'A'}{'score'} + $hash2{$scaffold}{$subseq}{'B'}{'score'},
		   );
		   next if $new_score < $current_score || $new_aln_length < $current_aln_length;
		   $hash2{$scaffold}{$subseq}{'A'} = $hash1{$subject}{$scaffold}{$subseq}{'A'};
		   $hash2{$scaffold}{$subseq}{'B'} = $hash1{$subject}{$scaffold}{$subseq}{'B'};		   
		}
	   }
	}
  }
  undef %hash1;

  # coz i'm too lazy to index but someone should
  my $longread_hash = &parse_fasta_to_hash($longreads);
  my $fill_in_hash = &parse_fasta_to_hash($fill_in);
  open (OUT,">$file.contigfill");
  open (OUTFSA,">$file.contigfill.fsa");
  foreach my $scaffold (keys %hash2){
	my $genome_seq = $fill_in_hash->{$scaffold};
	foreach my $subseq (sort {$b <=> $a} keys %{$hash2{$scaffold}}){
		if (
		    $hash2{$scaffold}{$subseq}{'A'} && $hash2{$scaffold}{$subseq}{'B'}
		    # now check if they are on same subject id; remember we already picked highest identity/score
		    # so may not have same subject ID
		    && $hash2{$scaffold}{$subseq}{'A'}{'subject'} eq $hash2{$scaffold}{$subseq}{'B'}{'subject'}
                    && $hash2{$scaffold}{$subseq}{'A'}{'orientation'} eq $hash2{$scaffold}{$subseq}{'B'}{'orientation'}
		    ){	
			my ($gap_replace_str,$orient,$subject) = ('',$hash2{$scaffold}{$subseq}{'A'}{'orientation'},$hash2{$scaffold}{$subseq}{'A'}{'subject'});

			if ($orient == 1){
				# Aend +1 -> Bstart  -1 but convert to 0-coord
				$gap_replace_str = substr( 
				  $longread_hash->{$subject},
				  $hash2{$scaffold}{$subseq}{'A'}{'s_real_end_1'} ,
				  $hash2{$scaffold}{$subseq}{'B'}{'s_real_start_1'} - $hash2{$scaffold}{$subseq}{'A'}{'s_real_end_1'} -1
				);
			}else{
				# Bend +1 ->Astart -1 but convert to 0-coord
				$gap_replace_str = substr( 
				  $longread_hash->{$subject},
				  $hash2{$scaffold}{$subseq}{'B'}{'s_real_end_1'} ,
				  $hash2{$scaffold}{$subseq}{'A'}{'s_real_start_1'} - $hash2{$scaffold}{$subseq}{'B'}{'s_real_end_1'} -1
				);
			}
			warn Dumper \$hash2{$scaffold}{$subseq};
#			die Dumper $gap_replace_str;
			print OUT $scaffold."\t".$hash2{$scaffold}{$subseq}{'A'}{'q_real_end_1'}
			  ."\t".$hash2{$scaffold}{$subseq}{'B'}{'q_real_start_1'}
			  ."\t".$hash2{$scaffold}{$subseq}{'A'}{'orig_gap_length'}
			  ."\t".length($gap_replace_str)
			  ."\n";
			# conv to 0-coord. replace gap and anything else we need to remove (not just orig_gap_length)

			# this will update the genome sequence so all the coordinates would change, that's why we replace from the last gap ;-)
			my $old_gap = substr($genome_seq,
					$hash2{$scaffold}{$subseq}{'A'}{'q_real_end_1'}-1,
					$hash2{$scaffold}{$subseq}{'B'}{'q_real_start_1'}-1 - $hash2{$scaffold}{$subseq}{'A'}{'q_real_end_1'}+1,
#					$hash2{$scaffold}{$subseq}{'A'}{'orig_gap_length'},
					$gap_replace_str );
warn Dumper ($old_gap,$gap_replace_str);

		}
	}
	print OUTFSA ">$scaffold\n".&wrap_text($genome_seq)."\n";
  }
 close OUT;
 close OUTFSA;
}


sub split_gap_fasta() {
 my $fasta = shift;
 my ($seq_counter,$subseq_counter) = (int(0),int(0));
 my $output = "$fasta.gapsplit.$gap_size_min.$gap_size_max.$capture_length";
 return $output if -s $output;
 open (OUT,">$output");
 my $orig_sep = $/;
 $/ = '>';
 open( IN, $fasta ) || die( "Cannot open $fasta : " . $! );
 while ( my $record = <IN> ) {
  chomp($record);
  next unless $record;
  my @lines = split( $orig_sep, $record );
  my $id    = shift(@lines);
  my $seq   = join( '', @lines );
  $seq =~ s/\s+//g;

  next unless $id && $seq;
  my $seq_length = length($seq);
  next unless $seq_length > 2 * $capture_length + $gap_size_min;
  $seq_counter++;
  while ($seq =~ /([^N]{$capture_length})(N{$gap_size_min,$gap_size_max})([^N]{$capture_length})/g ) {
   $subseq_counter++;
   my ($frag_left,$gap,$frag_right) =  ($1, $2, $3 );
   next unless $frag_left && $gap && $frag_right;
   #  5.6.0  The built-in variables @- and @+ hold the start and end positions, respectively, of the last successful match. 
   # is the offset of the start of the last successful match as  used in substring
   #$-[0] and $+[0] correspond to entire pattern, while $-[N] and $+[N] correspond to the $N ($1, $2) submatches)
   # NB: coordinates are from 0 so can use substr
	# so that gap : die Dumper (substr($seq,$-[2],$+[2]- $-[2]));
   my ($whole_start_pos,$whole_end_pos, 
	$left_start_pos, $left_end_pos,
	$right_start_pos, $right_end_pos,
	$gap_length
	) = ( $-[0], $+[0], 
	      $-[1], $+[1],
	      $-[3], $+[3],
	      ($+[2] - $-[2] ) );
   
   # so i don't forget: +1
   my ($whole_start_pos_1,$whole_end_pos_1) = ($whole_start_pos + 1, $whole_end_pos );
   my ($left_start_pos_1,$left_end_pos_1) = ($left_start_pos + 1, $left_end_pos );
   my ($right_start_pos_1,$right_end_pos_1) = ($right_start_pos + 1, $right_end_pos );

   # left
   print OUT ">$id:$subseq_counter".":A:$left_start_pos_1:$left_end_pos_1:$gap_length".$orig_sep
	.&wrap_text(substr($seq,$left_start_pos,$left_end_pos-$left_start_pos))
	.$orig_sep;
   # right
   print OUT ">$id:$subseq_counter".":B:$right_start_pos_1:$right_end_pos_1:$gap_length".$orig_sep
	.&wrap_text(substr($seq,$right_start_pos,$right_end_pos-$right_start_pos))
	.$orig_sep;
  }
 }
 close IN;
 close OUT;
 $/ = $orig_sep;
 print "Processed $seq_counter sequences into $subseq_counter pair of subsequences\n";
 return $output;
}

#############################################################
sub wrap_text() {
# there seems to be an error sometimes.
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 80 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/g;
 $string =~ s/\s+$//;
 return $string;
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}
###
sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}
