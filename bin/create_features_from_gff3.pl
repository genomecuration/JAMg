#!/usr/bin/env perl

=pod

=head1 NAME

 create_features_from_gff3.pl

=head1 USAGE

Provide a GFF3 file and a genome FASTA file to phase and create sequence feature (CDS,PEP,MRNA,GENE) information for
 
 Mandatory
  -gff|infile   :s GFF file to process
  -genome|fasta :s FASTA file of genome
 
 Optional
  -name               Use Transcript common name as the main ID in the output.
  -lettername         Use -R[two letters] notation for alternative transcript instead of .[digits]  
  -verbose            Print progress and debug info
  -skip_delete        Skip Status=Delete mRNAs
  -one_isoform        Process only one isoform per gene
  -rename             Rename the IDs with a new JAMg IDs, up to 32 characters (needed for hhblits). Don't use with -name
  -strip_name         Remove Name tag
  -change_source  :s  Change GFF Source to this value
  -simple             Don't add introns and splice sites in GFF
  -delete_ns      :i  If more than this many Ns, 1) remove UTRs from mRNA
  -genbank        :s  Prepare for submission to GenBank; give BioProject-assigned locus_tag_id 
  -go_file        :s  SQL query from JAMp database
  -fix_first_phase    Change phase of first CDS to 0 if it is not 0

NB: -name means that the common name has no spaces and is unique (will be checked). Useful for WebApollo

=head1 JAMp database

 To get the GO IDs you can run this query, e.g for dataset_1:

 SELECT ic.cds_uname,go_id,known_protein_id,significance,identity from dataset_1.inference_cds ic \
    JOIN known_proteins.go_assoc ga ON ic.known_protein_id=ga.uniprot_id  \
    ORDER BY significance ASC,identity DESC,character_length(known_protein_id),known_protein_id ASC

 Only first three columns are currently used. You can also define WHERE statements for significance or identity
 
=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Carp;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use GTF_utils;
use Nuc_translator;

$|=1;
our $SEE;
my $codon_table = 'universal';
my $minorf = 96;    #minimum orf size in bp
my $to_bed_exec = "$RealBin/../3rd_party/PASA/misc_utilities/gene_gff3_to_bed.pl";
my ($simple_gff3, $gfffile, $genome, $change_name,$lettername,$verbose, $fix_first_phase,
 $one_iso, $do_rename, $change_source, $strip_name, $split_single, $delete_ns, $bioproject_locus_id,$go_file,%dbxref_hash );
pod2usage $! unless &GetOptions(
            'one_isoform'      => \$one_iso,
            'gff|infile:s'     => \$gfffile,
            'genome|fasta:s'   => \$genome,
            'name|change_name' => \$change_name,
            'debug'            => \$SEE,
            'verbose'          => \$verbose,
            'lettername'       => \$lettername,
            'change_source:s'    => \$change_source,
      	    'rename'	    	   => \$do_rename,
            'strip_name'       => \$strip_name,
	    'simple' => \$simple_gff3,
	    'split_single' => \$split_single,
	    'delete_ns:i' => \$delete_ns,
	    'genbank:s' => \$bioproject_locus_id,
	    'go_file:s' =>\$go_file,
	    'codon_table:s' => \$codon_table,
	    'fix_first_phase' => \$fix_first_phase
);
$gfffile = shift if !$gfffile;
$genome  = shift if !$genome;
pod2usage unless $gfffile && $genome;

die "GFF file $gfffile does not exist\n"  unless -s $gfffile;
die "Fasta file $genome does not exist\n" unless -s $genome;

#global
Nuc_translator::use_specified_genetic_code ($codon_table);

if ($go_file){
	open (GO,$go_file) || die "Cannot find $go_file";
	while (my $ln=<GO>){
		chomp($ln);
		next if !$ln;
		my @data = split("\t",$ln);
		next unless $data[1] && $data[1]=~/^\d+$/;
		#hash hashes to prevent duplicates.
		$dbxref_hash{$data[0]}{'gos'}{'GO:'.$data[1]}++;
		$dbxref_hash{$data[0]}{'hits'}{'similar to AA sequence:UniProtKB:'.$data[2]}++ if $data[2];
	}
	close GO;
}
my $scaffold_seq_hashref = &read_fasta($genome);

my %unique_names_check;

my %params;
$params{unspliced_transcript} = 1;    # highlights introns

&gff3_process($gfffile);

#########################################################################
sub gff3_process() {
 my $gff3_file = shift;
 confess "Cannot find $gff3_file\n" unless -s $gff3_file;
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer;
$gene_obj_indexer =  new Gene_obj_indexer( { "create" => $index_file } );
 confess "Cannot index with $index_file\n" unless $gene_obj_indexer && -s $index_file;
 my $genome_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );

 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 open( GFF3, ">$gff3_file.gff3" );
 open( GFF3_ERR, ">$gff3_file.gff3.err" );
 open( GFF3_SINGLE, ">$gff3_file.gff3.single" ) if $split_single;
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );
 open( MRNA, ">$gff3_file.mRNA" );
 open( TRACK, ">$gff3_file.track"); 
 print "Indexing complete\n";
 my $gene_count = 1;
 if ($bioproject_locus_id){
	print GFF3 "##gff-version 3\n";
	print GFF3_SINGLE "##gff-version 3\n" if $split_single;
 }


 foreach my $genome_id (sort {
	($a =~ /(\d*)$/)[0] cmp ($b =~ /(\d*)$/)[0] || $a cmp $b
  } keys %$genome_id_to_gene_list_href ) {

  my $genome_seq = $scaffold_seq_hashref->{$genome_id};
  my $genome_seq_length = length($genome_seq);
  if ( !$genome_seq ) {
   warn "Cannot find genome sequence $genome_id\n";
   next;
  }

  my @gene_ids = sort {
    ($a =~ /(\d*)$/)[0] cmp ($b =~ /(\d*)$/)[0] || $a cmp $b
    } @{ $genome_id_to_gene_list_href->{$genome_id} };

  if ($bioproject_locus_id){
	if ($genome_seq=~/^N+/ || $genome_seq=~/N+$/){
		die "Sequence $genome_id begins or ends with a run of Ns (assembly gaps). This is not allowed for genome submissions, please edit the genome sequence\n";
	}elsif($genome_seq=~/([^ATCGN])/){
		die "Genome sequence sequence $genome_id has non-ATCGN characters (e.g. $1). Please run cleanup_iupac.pl or similar\n";
	}elsif ($genome_seq_length < 200){
		die "Sequence $genome_id is shorter than the GenBank minimum of 200. Please edit the genome sequence e.g. trim_fasta_all.pl -le 200 (or better -le 2000; I mean, are these even meaningful contigs?!)\n";
	}
	my $region_print = "##sequence-region\t$genome_id\t1\t$genome_seq_length\n";
	$region_print .= "$genome_id\tAssembly_$bioproject_locus_id\tregion\t1\t$genome_seq_length\t.\t+\t.\tID=$genome_id\n";
	print GFF3 $region_print;
	print GFF3_SINGLE $region_print if $split_single;
  }

  print "\nProcessing scaffold $genome_id\n" if $verbose;

 foreach my $gene_id (sort @gene_ids) {
    
   my (%preferences,$is_single_exon);
   $preferences{'sequence_ref'} = \$genome_seq;
   $preferences{'source'}  = $change_source if $change_source;

   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   next if ($gene_id=~/temp_model/  
	|| ($gene_obj_ref->{gene_name} && $gene_obj_ref->{gene_name}=~/temp_model/) 
	|| ($gene_obj_ref->{TU_feat_name} && $gene_obj_ref->{TU_feat_name}=~/temp_model/)
	);

   my ( $gene_lend, $gene_rend ) = sort { $a <=> $b } $gene_obj_ref->get_gene_span();
   if ($gene_rend > $genome_seq_length){
	warn "Skipping gene $gene_id as it expands past length of scaffold\n";
   }

   $gene_obj_ref->{TU_feat_name} = $gene_id if !$gene_obj_ref->{TU_feat_name};

   my $jamg_id_gene  = "JAMg_model_".$gene_count;
   my $old_gene_name = $gene_id;

   $gene_obj_ref->{genbank_submission} = $bioproject_locus_id if $bioproject_locus_id;
   $gene_obj_ref->create_all_sequence_types( \$genome_seq, %params );
   my $gene_seq = $gene_obj_ref->get_gene_sequence();

   next if !$gene_seq;
   $gene_seq =~ s/(\S{80})/$1\n/g;
   chomp $gene_seq;

   if ($do_rename){
     $gene_id = $jamg_id_gene;
     $gene_obj_ref->{gene_name} = $gene_id;
     $gene_obj_ref->{TU_feat_name} = $gene_id;
   }


   # check each isoform
   my (@isoforms,%peplength_hash);
   foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms() ){
	next unless $isoform;
	next if (
		($isoform->{Model_feat_name} && $isoform->{Model_feat_name} =~/temp_model/) 
		|| ( $isoform->{transcript_name} && $isoform->{transcript_name}=~/temp_model/)
	);
#warn Dumper (0,$isoform);

	my $do_recreate_seqs;
        my $transcript_main_id = $isoform->{Model_feat_name} || confess;
	my $error_name_text = "Transcript $transcript_main_id of $old_gene_name gene (new gene ID: $gene_id)";
	my $cds_seq = $isoform->get_CDS_sequence();
	my $mrna_seq = $isoform->get_cDNA_sequence();
	my $pep_seq = $isoform->get_protein_sequence();

	$isoform = &check_phasing($isoform,$error_name_text,\$genome_seq,1) unless $pep_seq!~/^M/; # see subroutine comments on why;


	if (!$pep_seq){
		warn "SKIP: $error_name_text has no protein translation.\n";
		next;
	}

	if ( $pep_seq!~/^M/ && $isoform->has_5prime_UTR()){
		warn "FIXING: $error_name_text has 5'UTR but does not start with Methionine. Removing 5' UTR\n";
		$isoform->trim_5UTR();
		$do_recreate_seqs++;
		$isoform->{is_5prime_partial} = 1;
	}

	if ($pep_seq!~/\*$/ && $isoform->has_3prime_UTR()){
		warn "FIXING: $error_name_text has 3'UTR but does not end with a stop codon. Removing 3' UTR\n";
		$isoform->trim_3UTR();
		$do_recreate_seqs++;
		$isoform->{is_3prime_partial} = 1;
	}
	if ($mrna_seq=~/^N/ && $isoform->has_5prime_UTR()){
		warn "FIXING: $error_name_text has 5'UTR endings in gaps. Removing 5' UTR\n";
		$isoform->trim_5UTR();
		$do_recreate_seqs++;
	}
	if ($mrna_seq=~/N$/ && $isoform->has_3prime_UTR()){
		warn "FIXING: $error_name_text has 3'UTR endings in gaps. Removing 3' UTR\n";
		$isoform->trim_3UTR();
		$do_recreate_seqs++;
	}

	if ($delete_ns){
		my $mNs = ($mrna_seq =~ tr/N//);
		if ($mNs && $mNs > $delete_ns){
			warn "UTR removed: $error_name_text has more than $delete_ns assembly gaps Ns ($mNs) in mRNA sequence.\n";
			$isoform = $isoform->trim_UTRs();
			$do_recreate_seqs++;
		}
   	}

	if ($do_recreate_seqs){
		$isoform->refine_gene_object();
		$isoform->create_all_sequence_types( \$genome_seq, %params );
		$mrna_seq = $isoform->get_cDNA_sequence();
		$cds_seq = $isoform->get_CDS_sequence();
		$gene_seq = $gene_obj_ref->get_gene_sequence();
		$pep_seq = $isoform->get_protein_sequence();
		$do_recreate_seqs = undef;
	}

	#always check for Ns in coding
	my $cNs = ($cds_seq =~ tr/N//);
	if ($cNs){
		#try to fix in beg or end; also necessary if bioproject genbank submission
		my ($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
		my @exons = $isoform->get_exons();
		#we need to do a while loop here because there could be concat exons with the same issue
		while ($cds_seq =~/^N/){
			my $adjust_length = ($cds_seq =~ tr/^N//);
			warn "FIXING: $error_name_text has a coding sequence starting with $adjust_length x non-standard nucleic code. Adjusting CDS to remove them\n";
			
			my $exon_coord5 = $exons[$first_cds_exon_number]->{end5};
			my $exon_coord3 = $exons[$first_cds_exon_number]->{end3};
			#set coords small to large
			my ( $ordered_exon_coord_start, $ordered_exon_coord_end ) = ( $exon_coord5 < $exon_coord3 ) ? ( $exon_coord5, $exon_coord3 ) : ( $exon_coord3, $exon_coord5 );
			my $exon_seq = substr( $genome_seq, ( $ordered_exon_coord_start - 1 ), ( $ordered_exon_coord_end - $ordered_exon_coord_start + 1 ) );

			if ($exon_seq=~/^N+$/){
				#then we have to delete it. we dont like empty elemnts in array so splice
				splice(@exons,$first_cds_exon_number,1);
			}else{
				#then we have to trim it - 5'
				if ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{'strand'} eq '+'){
					$exons[$first_cds_exon_number]->{end5} = $exon_coord5 + $adjust_length;
				}else{
					$exons[$first_cds_exon_number]->{end5} = $exon_coord5 - $adjust_length;# reversed coords...
				}
				$exons[$first_cds_exon_number]->{CDS_exon_obj}->{end5} = $exons[$first_cds_exon_number]->{end5};
				#5' trim means new phase
				$exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase} = $adjust_length % 3  ;
				#THIS wont work if partial: (length($cds_seq)-$adjust_length) % 3;
			}
			$isoform->{mRNA_exon_objs} = 0;
		   	$isoform->{mRNA_exon_objs} = \@exons;
			$isoform->refine_gene_object();
			$isoform->create_all_sequence_types( \$genome_seq, %params );
			($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
			$mrna_seq = $isoform->get_cDNA_sequence();
			$cds_seq = $isoform->get_CDS_sequence();
			$gene_seq = $gene_obj_ref->get_gene_sequence();
			$pep_seq = $isoform->get_protein_sequence();
		}
		#another while loop 
		while ($cds_seq =~/N$/){
			my $adjust_length = ($cds_seq =~ tr/N$//);
			warn "FIXING: $error_name_text has a coding sequence ending with $adjust_length x non-standard nucleic code. Adjusting CDS to remove them\n";
			my $exon_coord5 = $exons[$last_cds_exon_number]->{end5};
			my $exon_coord3 = $exons[$last_cds_exon_number]->{end3};
			#set coords small to large
			my( $ordered_exon_coord_start, $ordered_exon_coord_end ) = ( $exon_coord5 < $exon_coord3 ) ? ( $exon_coord5, $exon_coord3 ) : ( $exon_coord3, $exon_coord5 );
			my $exon_seq = substr( $genome_seq, ( $ordered_exon_coord_start - 1 ), ( $ordered_exon_coord_end - $ordered_exon_coord_start + 1 ) );
			if ($exon_seq=~/^N+$/){
				#then we have to delete it. we dont like empty elemnts in array so splice
				splice(@exons,$last_cds_exon_number,1);
			}else{
				#then we have to trim it 3'
				if ($exons[$last_cds_exon_number]->{CDS_exon_obj}->{'strand'} eq '+'){
					$exons[$last_cds_exon_number]->{end3} = $exon_coord3 - $adjust_length;
				}else{
					$exons[$last_cds_exon_number]->{end3} = $exon_coord3 + $adjust_length;# reversed coords...
				}
				$exons[$last_cds_exon_number]->{CDS_exon_obj}->{end3} = $exons[$last_cds_exon_number]->{end3};
			}
			$isoform->{mRNA_exon_objs} = 0;
		   	$isoform->{mRNA_exon_objs} = \@exons;
			$isoform->refine_gene_object(); 
  			$isoform->create_all_sequence_types(\$genome_seq, %params);
			($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
			$mrna_seq = $isoform->get_cDNA_sequence();
			$cds_seq = $isoform->get_CDS_sequence();
			$gene_seq = $gene_obj_ref->get_gene_sequence();
			$pep_seq = $isoform->get_protein_sequence();
		}
		$cNs = ($cds_seq =~ tr/N//);
	}
#warn Dumper (1,$isoform);
	if ($cNs){
		if ( ($cNs / length($cds_seq) >= 0.5) ){
			warn "SKIP: $error_name_text has assembly gaps Ns ($cNs) in coding sequence that are more than 50%. Will print as error\n";
			print GFF3_ERR $isoform->to_GFF3_format(%preferences) . "\n";
			next;
		}
		warn "WARNING: $error_name_text has assembly gaps Ns ($cNs) in coding sequence.\n";
		if (!$bioproject_locus_id){
			print GFF3_ERR $isoform->to_GFF3_format(%preferences) . "\n";
			next;
		}
		$isoform->{internal_gap} = 1;
	}
#warn Dumper (2,$isoform);
	if ($fix_first_phase){
		 # ah, we do set the first phase to > 0 if we are next to a gap/splice site and want to maintain the frame (later)
 		my ($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	        my @exons = $isoform->get_exons();
		# finally if phase of first cds >0 then we need to chop and reset
		my $first_cds = $exons[$first_cds_exon_number]->get_CDS_obj();
	        my $first_phase = $first_cds->{phase};
		  if ( ref $first_cds && $first_phase > 0 ) {
			warn "FIXING: $error_name_text has phase of first CDS > 0 ($first_phase), setting to 0 and resetting co-ordinates by $first_phase. This may change again later though, stay tuned\n";
	                if ($first_cds->{strand} eq '+'){
	                        # might have utr so only if no utr
	                        $exons[$first_cds_exon_number]->{end5}+=$first_phase if $exons[$first_cds_exon_number]->{end5} == $first_cds->{end5};
	                        $first_cds->{end5}+=$first_phase;
	                }else{
	                        $exons[$first_cds_exon_number]->{end5}-=$first_phase if $exons[$first_cds_exon_number]->{end5} == $first_cds->{end5};
	                        $first_cds->{end5}-=$first_phase;
	                }
			$first_cds->{phase} = int(0);
		        $first_phase = int(0);
		   }
		$isoform->{mRNA_exon_objs} = 0;
	   	$isoform->{mRNA_exon_objs} = \@exons;
		$isoform->refine_gene_object();
		$isoform->create_all_sequence_types( \$genome_seq, %params );
		$mrna_seq = $isoform->get_cDNA_sequence();
		$cds_seq = $isoform->get_CDS_sequence();
		$gene_seq = $gene_obj_ref->get_gene_sequence();
		$pep_seq = $isoform->get_protein_sequence();
	}


	# sometimes breaks if we allow phase check with partial 5' genes but we have to do this because pasa rephases evm the wrong way
	$isoform = &check_phasing($isoform,$error_name_text,\$genome_seq,1) unless $pep_seq!~/^M/;
#warn Dumper (3,$isoform);

	$mrna_seq = $isoform->get_cDNA_sequence();
	$cds_seq = $isoform->get_CDS_sequence();
	$gene_seq = $gene_obj_ref->get_gene_sequence();
	$pep_seq = $isoform->get_protein_sequence();

	my $strand = $isoform->{strand};
	my @gene_span = $isoform->get_gene_span();

	# EXTEND IF PARTIAL + 1-3 NUCS FROM END - no phase check after. 
	if ($gene_span[0] < 4 && $gene_span[0] != 1){
		my @exons = $isoform->get_exons();
		if ( $pep_seq!~/^M/ && $strand eq '+' ){
			#expand gene, first exon, its cds if there
			$exons[0]{'end5'} = 1;
			if ( $exons[0]{'CDS_exon_obj'} && $exons[0]{'CDS_exon_obj'}{'end5'}){
				$exons[0]{'CDS_exon_obj'}{'end5'} = $exons[0]{'end5'};
			}
		}elsif($pep_seq!~/\*$/ && $strand eq '-'){
			#expand gene, last exon, its cds if there
			$exons[-1]{'end5'} = 1;
			if ( $exons[-1]{'CDS_exon_obj'} && $exons[-1]{'CDS_exon_obj'}{'end5'}){
				$exons[-1]{'CDS_exon_obj'}{'end5'} = $exons[-1]{'end5'};
			}
		}else{
			next;
		}
		warn "FIX for GenBank: $error_name_text has partial CDS very near the 5' end of the genome sequence\n";
		$isoform->{mRNA_exon_objs} = 0;
	   	$isoform->{mRNA_exon_objs} = \@exons;
		$isoform->refine_gene_object(); 
		$isoform->create_all_sequence_types( \$genome_seq, %params );
		$mrna_seq = $isoform->get_cDNA_sequence();
		$cds_seq = $isoform->get_CDS_sequence();
		$gene_seq = $gene_obj_ref->get_gene_sequence();
		$pep_seq = $isoform->get_protein_sequence();
		@gene_span = $isoform->get_gene_span();
	}
#warn Dumper (4,$isoform);
	if ($gene_span[1] > $genome_seq_length - 4 && $gene_span[1] != $genome_seq_length){
		my @exons = $isoform->get_exons();
		if ( $pep_seq!~/^M/ && $strand eq '-' ){
			#expand gene, first exon, its cds if there
			$exons[0]{'end3'} = $genome_seq_length;
			if ( $exons[0]{'CDS_exon_obj'} && $exons[0]{'CDS_exon_obj'}{'end3'}){
				$exons[0]{'CDS_exon_obj'}{'end3'} = $exons[0]{'end3'};
			}
		}elsif($pep_seq!~/\*$/ && $strand eq '+'){
			#expand gene to $genome_seq_length, last exon, its cds if there
			$exons[-1]{'end3'} = $genome_seq_length;
			if ( $exons[-1]{'CDS_exon_obj'} && $exons[-1]{'CDS_exon_obj'}{'end3'}){
				$exons[-1]{'CDS_exon_obj'}{'end3'} = $exons[-1]{'end3'};
			}
		}else{
			next;
		}
		warn "FIX for GenBank: $error_name_text has partial CDS very near the 3' end of the genome sequence\n";
		$isoform->{mRNA_exon_objs} = 0;
	   	$isoform->{mRNA_exon_objs} = \@exons;
		$isoform->refine_gene_object(); 
		$isoform->create_all_sequence_types(\$genome_seq, %params);
		$mrna_seq = $isoform->get_cDNA_sequence();
		$cds_seq = $isoform->get_CDS_sequence();
		$gene_seq = $gene_obj_ref->get_gene_sequence();
		$pep_seq = $isoform->get_protein_sequence();
		@gene_span = $isoform->get_gene_span();
	}
#warn Dumper (5,$isoform);

	# the following section MUST not be followed by a phase correction to 0

	if (length($cds_seq) % 3 != 0 && $pep_seq!~/\*$/){
		warn "FIXING: $error_name_text has no stop codon, a CDS that is (length ".length($cds_seq)." whose % 3 != 0). Scanning for a stop codon\n";
		$isoform = &find_next_stop($isoform,$error_name_text,\$genome_seq,1);
	}elsif($pep_seq!~/\*$/){
		warn "FIXING: $error_name_text has no stop codon. Scanning for a stop codon\n";
		$isoform = &find_next_stop($isoform,$error_name_text,\$genome_seq);
	}
#warn Dumper (6,$isoform);
	if (length($cds_seq) % 3 != 0 && $pep_seq!~/^M/){
		warn "FIXING: $error_name_text has no start codon, a CDS that is (length ".length($cds_seq)." whose % 3 != 0). Scanning for a start codon\n";
		$isoform = &find_next_start($isoform,$error_name_text,\$genome_seq,1);
	}elsif ($pep_seq!~/^M/){
		warn "FIXING: $error_name_text has no start codon. Scanning for a start codon\n";
		$isoform = &find_next_start($isoform,$error_name_text,\$genome_seq);
	}
	# Do not do any more phase correction to 0
#warn Dumper (7,$isoform);
	$mrna_seq = $isoform->get_cDNA_sequence();
	$cds_seq = $isoform->get_CDS_sequence();
	$gene_seq = $gene_obj_ref->get_gene_sequence();
	$pep_seq = $isoform->get_protein_sequence();
	@gene_span = $isoform->get_gene_span();

	$peplength_hash{$transcript_main_id} = length($pep_seq);
    	if ( length($pep_seq) < int($minorf/3) ) {
		warn "SKIP: $error_name_text coding sequence is shorter than $minorf base pairs.\n";
		next;
	}
	if ($pep_seq=~/\*\S/){
	        warn "WARNING: $error_name_text has stop codon inside ORF. If a GenBank submission, then then tagging as pseudo (broken assembly)\n";
		$isoform->{internal_stop} = 1;
		if (!$bioproject_locus_id){
			print GFF3_ERR $isoform->to_GFF3_format(%preferences) . "\n";
			next;
		}
	}
#warn Dumper (8,$isoform);


	# if we have survived this far:
	if ($one_iso){
		#pick longest
		$isoforms[0] = $isoform if !$isoforms[0] || $peplength_hash{$transcript_main_id} > $peplength_hash{$isoforms[0]->{Model_feat_name}};
		
	}else{
		push(@isoforms,$isoform); # all data
	}
   }


   next unless $isoforms[0];

   if ($one_iso){
	$isoforms[0]->delete_isoforms();
   }

   if (scalar(@isoforms)>1){
	my @iso_check;
	# now check if any isoforms are duplicates of each other. Duplicates are defined as same exon coords
	for (my $i=0;$i<scalar ( @isoforms );$i++){
		next if !$isoforms[$i];
		my @exonsI = $isoforms[$i]->get_exons();
		for (my $k=1;$k<scalar ( @isoforms ) ;$k++){
			next if $i==$k;
			if (!$isoforms[$k]){
				next;
			}
			my $safe_flag;
			my @exonsK = $isoforms[$k]->get_exons();
			next if scalar(@exonsI) != scalar(@exonsK);
			for (my $x=0;$x<scalar(@exonsI);$x++){
				my @coordI = $exonsI[$x]->get_coords();
				my @coordK = $exonsK[$x]->get_coords();
				if ($coordI[0] != $coordK[0] || $coordI[1] != $coordK[1]){
					$safe_flag++;
					last;
				}
			}
			if (!$safe_flag){
				my $transcript_main_idI = $isoforms[$i]->{Model_feat_name};
				my $transcript_main_idK = $isoforms[$k]->{Model_feat_name};
				warn "Isoforms $transcript_main_idI and $transcript_main_idK of $old_gene_name ($gene_id new gene name) are identical. Deleting second one.\n";
				for (my $x=0;$x<scalar(@{$isoforms[$i]{additional_isoforms}});$x++){
					delete ($isoforms[$i]{additional_isoforms}[$x]) if ($isoforms[$i]{additional_isoforms}[$x] && $isoforms[$i]{additional_isoforms}[$x]->{Model_feat_name} eq $transcript_main_idK);
				}
				delete($isoforms[$k]); # so we don't process it again
			}
		}
		push(@iso_check,$isoforms[$i]);
	}
	@isoforms = @iso_check;
   }

   # verify the gene span is not larger than any isoform span
   my @gene_span = $gene_obj_ref->{gene_span};
   my @new_gene_span;
   

   foreach my $isoform ( @isoforms )   {

	my @exons = $isoform->get_exons();
	foreach my $e (@exons){
		my @coords = $e->get_coords();
		if ($isoform->{strand} eq '+'){
			$new_gene_span[0] = $coords[0] if (!$new_gene_span[0] || $new_gene_span[0] > $coords[0]);
			$new_gene_span[1] = $coords[1] if (!$new_gene_span[1] || $new_gene_span[1] < $coords[1]);
		}else{
			$new_gene_span[0] = $coords[0] if (!$new_gene_span[0] || $new_gene_span[0] < $coords[0]);
			$new_gene_span[1] = $coords[1] if (!$new_gene_span[1] || $new_gene_span[1] > $coords[1]);
		}
	}
  }
  $gene_obj_ref->{gene_span}   = [ $new_gene_span[0], $new_gene_span[1] ];

   if ($do_rename){
     print TRACK "GENE\t$gene_id\t$jamg_id_gene\n";
   }
   print GENE ">$gene_id type:gene\n$gene_seq\n";
   print "\rprocessing gene $old_gene_name ($gene_id new gene name)"
     . "                                                  "
     if $verbose;

   my $isoform_count;
   foreach my $isoform ( @isoforms )   {
    my $isoform_id  = $isoform->{Model_feat_name};

    next unless $isoform->has_CDS() || !$isoform->get_CDS_span();
    my @model_span  = $isoform->get_CDS_span();
    next if ( abs( $model_span[0] - $model_span[1] ) < 3 );
    
   # there is a real issue here. Most times EVM and PASA are ok with phasing
   # but sometimes there are issues.
   # and these issues happen both ways (with and without rephasing). 
   # so i decided to do rephasing always 
   # unfortunately this has to happen across the entire gene
   # and there are 2-3 genes (out of 20000) that have issues.
   # The second issue was start CDSs that had phase !=0 
   # this is now solved within the library via set_CDS_phases (basically first CDS is reset to have phase = 0 (frame = 1)
	# JAN19; not solved.... it seems it was only for gtf
   if ($isoform->{'num_exons'} == 1){
	$is_single_exon++;
   }
    print "\rprocessing gene $old_gene_name ($gene_id new gene name) isoform $isoform_id                                       " if $verbose;
    $isoform_count++;
    my $common_name = $isoform->{transcript_name} || $isoform->{com_name};
    my $description = '';
    my $alt_name    = '';
    my $transcript_main_id     = $isoform_id;
    
    if ( $common_name && $change_name ) {
     if ( $common_name =~ /\s/ ) {
      $description = $common_name;
      $description =~ s/^\s*(\S+)\s*//;
      $description = $description;
      $common_name = $1 || die;
     }

     if ( $common_name =~ /\.\d+$/ && !$lettername){ 
      $transcript_main_id = $common_name;
     }elsif ($lettername && $common_name =~ /-R[A-Z]+$/ ) {
      $transcript_main_id = $common_name;
     }elsif ($lettername){
       $transcript_main_id = $common_name . '-RA';
     }else{
       $transcript_main_id = $common_name . '.1';
     }
     $alt_name = "($isoform_id) ";


     if ( $unique_names_check{$common_name} ) {
      if ($lettername && $common_name =~ /(-R[A-Z]+)$/ ) {
       warn "Common name $common_name ends in transcript notation ($1) but it is not unique!\n";
      }
      elsif (!$lettername && $common_name =~ /(\.\d+)$/) {
       warn "Common name $common_name ends in transcript notation ($1) but it is not unique!\n";
      }
      if ($lettername){
       my $letter = 'B';
       for (my $i=1;$i<$unique_names_check{$common_name};$i++){
        $letter++;
       }
       $transcript_main_id  = $common_name . '-R'.$letter;
      }else{
       $transcript_main_id  = $common_name . '.'.($unique_names_check{$common_name}+1);
      }
      $alt_name = "($isoform_id) ";
     }
     $unique_names_check{$common_name}++;
     print TRACK "TRANSCRIPT\t$isoform_id\t$transcript_main_id\n";
     # set description as note and update name
     $isoform->{transcript_name} = $transcript_main_id;
     $isoform->{com_name} = $transcript_main_id;
     $isoform->{pub_comment} = $description if $description;
    } # end change name
    elsif($do_rename){
      print TRACK "TRANSCRIPT\t$isoform_id\t$jamg_id_gene.$isoform_count\n";
      $transcript_main_id = $jamg_id_gene.'.'.$isoform_count;
      $isoform->{Model_feat_name} = $transcript_main_id;
    }

    if ($strip_name){
      $isoform->{transcript_name} = $transcript_main_id;
      $isoform->{com_name} = $transcript_main_id;
    }

    # get sequences
    my $pep_seq = $isoform->get_protein_sequence();
    my $cds_seq = $isoform->get_CDS_sequence();
    my $mrna_seq = $isoform->get_cDNA_sequence();
    # proteins
    $pep_seq =~ s/(\S{80})/$1\n/g;
    chomp $pep_seq if $pep_seq;
    print PEP ">$transcript_main_id ".$alt_name."type:polypeptide  gene:$gene_id$description\n".uc($pep_seq)."\n";

    # CDS
    $cds_seq =~ s/(\S{80})/$1\n/g if $cds_seq;
    chomp $cds_seq if $cds_seq;
    print CDS ">$transcript_main_id ".$alt_name."type:CDS  gene:$gene_id$description\n".uc($cds_seq)."\n";

    # mRNA (all exons)
    $mrna_seq =~ s/(\S{80})/$1\n/g;
    chomp $mrna_seq;
    print MRNA ">$transcript_main_id ".$alt_name."type:mRNA gene:$gene_id$description\n".uc($mrna_seq)."\n";

      if ($dbxref_hash{$transcript_main_id}){
	@{$isoform->{Dbxref}} = keys %{$dbxref_hash{$transcript_main_id}{'gos'}};
	@{$isoform->{inference}} = keys %{$dbxref_hash{$transcript_main_id}{'hits'}};
      }

   }


   # GFF3
   if ($simple_gff3){
	if ($is_single_exon && $split_single){
	        print GFF3_SINGLE $gene_obj_ref->to_GFF3_format(%preferences) . "\n";
	}else{
	        print GFF3 $gene_obj_ref->to_GFF3_format(%preferences) . "\n";
	}
   }else{
	if ($is_single_exon && $split_single){
	   	print GFF3_SINGLE $gene_obj_ref->to_GFF3_format_extended(%preferences) . "\n";
	}else{
	   	print GFF3 $gene_obj_ref->to_GFF3_format_extended(%preferences) . "\n";
	}
   }
   $gene_count++;
  }
 }
 close GFF3;
 close GFF3_ERR;
 close GFF3_SINGLE;
 unlink("$gff3_file.gff3.err") if !-s "$gff3_file.gff3.err" || -s "$gff3_file.gff3.err" < 100;
 close PEP;
 close CDS;
 close MRNA;
 close GENE;
 close TRACK;

 #&sort_gff3("$gff3_file.gff3");
 unlink("$gff3_file.track") unless -s "$gff3_file.track";

 unlink("$gff3_file.gff3.single") if !-s "$gff3_file.gff3.single" || -s "$gff3_file.gff3.single" < 100;
 #&sort_gff3("$gff3_file.gff3.single") if -s "$gff3_file.gff3.single";

 if ($to_bed_exec && -x $to_bed_exec){
	system("$to_bed_exec $gff3_file.gff3".' | sed -r \'~s/;\S+//\' | sed \'~s/ID=//\' '. " > $gff3_file.bed");
 }

 return;
}

sub read_fasta() {
 my $fasta = shift;
 my %hash;
 my $orig_sep = $/;
 $/ = '>';
 open( IN, $fasta ) || confess( "Cannot open $fasta : " . $! );
 while ( my $record = <IN> ) {
  chomp($record);
  next unless $record;
  my @lines = split( "\n", $record );
  my $id    = shift(@lines);
  my $seq   = join( '', @lines );
  $seq =~ s/\s+//g;
  if ( $id && $seq && $id =~ /^(\S+)/ ) {
   $hash{$1} = $seq;
  }
 }
 close IN;
 $/ = $orig_sep;
 return \%hash;
}

sub sort_gff3() {

 # find out why two empty lines still cause sort to crash
 my $gff       = shift;
 my $delimiter = shift;
 $delimiter = &get_gff_delimiter($gff) if !$delimiter;
 my $orig_sep = $/;
 open( GFF, $gff ) || confess "Cannot open $gff $!";
 open( GFF3, ">$gff.sorted" );
 $/ = $delimiter;
 my @records;

 while ( my $rec = <GFF> ) {
  chomp($rec);
  next if !$rec || $rec =~ /^\s*$/;
  push( @records, $rec );
 }
 close GFF;

 foreach my $rec (
  sort {
   my @first  = split( "\n", $a );
   my @second = split( "\n", $b );

   my @split1 = split( "\t", $first[0] );
   my @split2 = split( "\t", $second[0] );
   return $split1[0] cmp $split2[0] || $split1[3] <=> $split2[3];
  } @records
   )
 {
  print GFF3 $rec . $/;
 }

 close GFF3;
 $/ = $orig_sep;
 rename( "$gff.sorted", $gff );
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

##########################################################3
sub fix_single_exon_phase(){  # not used
	my $obj = shift;
	my $strand = $obj->get_orientation();
	 foreach my $exon ( $obj->get_exons() ) {
	  if ( my $cds = $exon->get_CDS_exon_obj() ) {
		my $phase = $cds->get_phase();
		next if !$phase || $phase == 0;
		#delete this many from beginning
		if ($strand eq '+'){
			$exon->{end5}+= $phase;
			$cds->{end5}+= $phase;
			$cds->{phase}= 0;			
		}else{
			# i think
			$exon->{end5}-= $phase;
			$cds->{end5}-= $phase;
			$cds->{phase}= 0;	
		}

	  }
	 }

	return $obj;
}
#######################################
sub find_next_stop(){
	my ($isoform,$error_name_text,$genome_seq_ref,$change_phase) = @_;
	my @stops = &Nuc_translator::get_stop_codons();
	my %stop_codons; 
	foreach my $s (@stops){$stop_codons{$s}++;}

	my $strand = $isoform->{strand};
	my @exons = $isoform->get_exons();
	my ($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	my $first_phase = $exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase};
	my $cds_seq = $isoform->get_CDS_sequence();
	my $cds_seq_length = length($cds_seq);
	my $cds_modulo_gtf_phase = $cds_seq_length % 3;
	# fun fun fun: GTF (and our libraries) work with a different definition
	# of phase, where it is which codon base the first base corresponds to
	# 0,1,2 for first, second, third base. 
	if ($cds_modulo_gtf_phase == 2){$cds_modulo_gtf_phase=1;}
	elsif($cds_modulo_gtf_phase == 1){$cds_modulo_gtf_phase=2;}

	my $phase_to_check = $change_phase ? $first_phase : $cds_modulo_gtf_phase;

	my ($c1,$c2) = $exons[$last_cds_exon_number]->get_coords();
	my ($next_two,$cds_last_codon);
	if ($strand eq '-' ){
		my $t = $c1;
		$c1 = $c2;
		$c2 = $t;
		$cds_last_codon = &reverse_complement(substr( $$genome_seq_ref,( $c1 - 1 ),3));
		$next_two = &reverse_complement(substr($$genome_seq_ref,$c1-3,2));
        }else{
		$cds_last_codon = substr( $$genome_seq_ref,( $c2 -3 ),3);
		$next_two = substr($$genome_seq_ref,$c2,2);
	}
	return $isoform if $next_two eq 'GT' || $next_two eq 'GC' || $next_two =~/^N/;
	return $isoform if $stop_codons{$cds_last_codon};

	#we need to make sure there is no 3'UTR
	$isoform->trim_3UTR();
	$isoform->refine_gene_object(); 
	$isoform->create_all_sequence_types( $genome_seq_ref, %params );
	#have to redo this as we removed UTR
	@exons = $isoform->get_exons();
	($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	($c1,$c2) = $exons[$last_cds_exon_number]->get_coords();
	if ($strand eq '-' ){
		my $t = $c1;
		$c1 = $c2;
		$c2 = $t;
		$cds_last_codon = &reverse_complement(substr( $$genome_seq_ref,( $c1 - 1 ),3));
		$next_two = &reverse_complement(substr($$genome_seq_ref,$c1-3,2));
        }else{
		$cds_last_codon = substr( $$genome_seq_ref,( $c2 -3 ),3);
		$next_two = substr($$genome_seq_ref,$c2,2);
	}

	#if the phase is wrong, carry on search || if phase is correct but you
	# still haven't found a stop codon, carry on search
	# last will force the exit on intron splice site, gap or end of search
	# natural exit if phase is correct and stop codon found
	while ($cds_modulo_gtf_phase != $phase_to_check || !$stop_codons{$cds_last_codon} ){
		$cds_seq_length++;
		if ($strand eq '-') {
			last if $c1==1;
			$c1--;	
			$cds_last_codon = &reverse_complement(substr( $$genome_seq_ref,( $c1 - 1 ),3));
			$next_two = &reverse_complement(substr($$genome_seq_ref,$c1-3,2));
		}else{
			last if $c2==length($$genome_seq_ref);
			$c2++;
			$cds_last_codon = substr( $$genome_seq_ref,( $c2 -3 ),3);
			$next_two = substr($$genome_seq_ref,$c2,2);
		}
		last if $next_two =~/^N/;
		if ($next_two eq 'GT' || $next_two eq 'GC'){
			last;
		}
		$cds_modulo_gtf_phase = $cds_seq_length % 3;
		if ($cds_modulo_gtf_phase == 2){$cds_modulo_gtf_phase=1;}
		elsif($cds_modulo_gtf_phase == 1){$cds_modulo_gtf_phase=2;}
	}
	if ($exons[$last_cds_exon_number]->{CDS_exon_obj}->{'strand'} eq '+'){
		$exons[$last_cds_exon_number]{end3} = $c2;
		$exons[$last_cds_exon_number]->{CDS_exon_obj}->{end3} = $c2;
	}else{
		$exons[$last_cds_exon_number]{end3} = $c1;
		$exons[$last_cds_exon_number]->{CDS_exon_obj}->{end3} = $c1;
	}
	$isoform->{mRNA_exon_objs} = 0;
   	$isoform->{mRNA_exon_objs} = \@exons;
	$isoform->refine_gene_object(); 
	$isoform->create_all_sequence_types( $genome_seq_ref, %params );
	return $isoform;
}

sub find_next_start(){
	my ($isoform,$error_name_text,$genome_seq_ref,$change_phase) = @_;
	my $cds_seq = $isoform->get_CDS_sequence();
	my $cds_seq_length = length($cds_seq);
	my $cds_modulo_gtf_phase = $cds_seq_length % 3;
#	# fun fun fun: GTF (and our libraries) work with a different definition
#	# of phase, where it is which codon base the first base corresponds to
#	# 0,1,2 for first, second, third base. 
	if ($cds_modulo_gtf_phase == 2){$cds_modulo_gtf_phase=1;}
	elsif($cds_modulo_gtf_phase == 1){$cds_modulo_gtf_phase=2;}
	my @stops = &Nuc_translator::get_stop_codons();
	my %stop_codons; 
	foreach my $s (@stops){$stop_codons{$s}++;}

	my @exons = $isoform->get_exons();
	my $strand = $isoform->{strand};
	my ($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	my $first_phase = $exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase};

	my ($c1,$c2) = $exons[$first_cds_exon_number]->get_coords();
	my ($next_two,$next_three,$cds_first_codon);
	if ($strand eq '-' ){
		my $t = $c1;
		$c1 = $c2;
		$c2 = $t;
		$cds_first_codon = &reverse_complement(substr( $$genome_seq_ref, $c2 -3,3));
		$next_two = &reverse_complement(substr($$genome_seq_ref,$c2,2));
		$next_three = &reverse_complement(substr($$genome_seq_ref,$c2,3));
        }else{
		$cds_first_codon = substr($$genome_seq_ref,$c1-1,3);
		$next_two = substr($$genome_seq_ref,$c1-1-2,2);
		$next_three = substr($$genome_seq_ref,$c1-1-3,3);
	}

        my $pep_seq = $isoform->get_protein_sequence();
	my $phase_to_check = $change_phase ? $first_phase : $cds_modulo_gtf_phase;
	my $is_3_partial = 1 if ($pep_seq!~/\*$/);
	# we have this issue: cds_seq_length % 3 cannot be used to define 
	# the right phase when we have a 3' partial gene. in that case cds_module must be used
	# this is not necessarily right since we don't know the true ORF/phase
	# but at least it prevents from the peptide sequence from changing.
	# in reality however, this search may not be the best, sometimes 
	# we should search for a start codon upstream but a splice site downstream
	$phase_to_check = $is_3_partial ? $cds_modulo_gtf_phase : $phase_to_check;	

	return $isoform if $next_two eq 'AG' || $next_two =~/N$/ || ($cds_seq_length % 3 == $phase_to_check && $stop_codons{$next_three});
	#we need to make sure there is no 5'UTR
	$isoform->trim_5UTR();
	$isoform->refine_gene_object(); 
	$isoform->create_all_sequence_types( $genome_seq_ref, %params );
	warn "WARNING: $error_name_text has an internal stop codon\n" if $pep_seq=~/\*[A-Z]/;

	#have to redo this as we removed UTR
	@exons = $isoform->get_exons();
	($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	($c1,$c2) = $exons[$first_cds_exon_number]->get_coords();
	if ($strand eq '-' ){
		my $t = $c1;
		$c1 = $c2;
		$c2 = $t;
		$cds_first_codon = &reverse_complement(substr( $$genome_seq_ref, $c2 -3,3));
		$next_two = &reverse_complement(substr($$genome_seq_ref,$c2,2));
		$next_three = &reverse_complement(substr($$genome_seq_ref,$c2,3));
        }else{
		$cds_first_codon = substr($$genome_seq_ref,$c1-1,3);
		$next_two = substr($$genome_seq_ref,$c1-1-2,2);
		$next_three = substr($$genome_seq_ref,$c1-1-3,3);
	}

	my $phase_adjustment = $first_phase;
	while ($cds_modulo_gtf_phase != $phase_to_check || $cds_first_codon ne 'ATG' ){
		$cds_seq_length++;
		if ($strand eq '-') {
			last if $c2==length($$genome_seq_ref);
			$c2++;
			$cds_first_codon = &reverse_complement(substr( $$genome_seq_ref, $c2 -3,3));
			$next_two = &reverse_complement(substr($$genome_seq_ref,$c2,2));
			$next_three = &reverse_complement(substr($$genome_seq_ref,$c2,3));
		}else{
			last if $c1==1;
			$c1--;
			$cds_first_codon = substr($$genome_seq_ref,$c1-1,3);
			$next_two = substr($$genome_seq_ref,$c1-1-2,2);
			$next_three = substr($$genome_seq_ref,$c1-1-3,3);
		}
		$phase_adjustment++;
		my $phase = $phase_adjustment % 3;
		if ($phase == 2){$phase=1;}
		elsif($phase == 1){$phase=2;}
		$cds_modulo_gtf_phase = $cds_seq_length % 3;
		if ($cds_modulo_gtf_phase == 2){$cds_modulo_gtf_phase=1;}
		elsif($cds_modulo_gtf_phase == 1){$cds_modulo_gtf_phase=2;}
		
		last if $cds_modulo_gtf_phase == $phase_to_check && $stop_codons{$next_three};
		# phase change needed? often we have the partial 5' being at a splice site
		# but originally erroneously marked by a UTR region
		# this does reset the first CDS to a new number
#warn Dumper ($phase_adjustment,$cds_first_codon,$phase);
		if ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase} != $phase && ($next_two eq 'AG' || $next_two =~/N$/)){
			warn "FIXING: $error_name_text has a new phase for the first CDS because it is adjacent to a splice site or gap"
#				."(".$exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase}
#				." vs "
#				.$phase .")"
				."\n";
			$exons[$first_cds_exon_number]->{CDS_exon_obj}->{phase} = $phase;
			last;
		}
	}

	if ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{'strand'} eq '+'){
		$exons[$first_cds_exon_number]{end5} = $c1;
		$exons[$first_cds_exon_number]->{CDS_exon_obj}->{end5} = $c1;
	}else{
		$exons[$first_cds_exon_number]{end5} = $c2;
		$exons[$first_cds_exon_number]->{CDS_exon_obj}->{end5} = $c2;
	}
	# ok, rewrite isoform object
	$isoform->{mRNA_exon_objs} = 0;
   	$isoform->{mRNA_exon_objs} = \@exons;
	$isoform->refine_gene_object(); 
	$isoform->create_all_sequence_types( $genome_seq_ref, %params );
        $pep_seq = $isoform->get_protein_sequence();
	warn "MANUAL FIX: $error_name_text has an internal stop codon\n" if $pep_seq=~/\*[A-Z]/;
	return $isoform;
}

sub get_cds_ends(){
	my $isoform = shift;

	my @exons = $isoform->get_exons();
	my ($first_cds,$first_exon,$last_cds,$last_exon,$first_cds_exon_number,$last_cds_exon_number);
	for (my $i=0;$i<scalar(@exons);$i++){
		$first_exon = $exons[$i];
		$first_cds = $first_exon->get_CDS_obj();
		$first_cds_exon_number = $i;
		last if $first_cds;
	}
	for (my $i=scalar(@exons)-1;$i>=0;$i--){
		$last_exon = $exons[$i];
		$last_cds = $last_exon->get_CDS_obj();
		$last_cds_exon_number = $i;
		last if $last_cds;
	}

	return ($first_cds_exon_number,$last_cds_exon_number);

}

sub check_phasing(){
	my ($isoform,$error_name_text,$genome_seq_ref,$do_phasing) = @_;
	return if !$isoform;
	#reset the phases based on the longest CDS within the ORF. doesn't change much but when it does, it 
	# completely destroys the CDS... this happens when the first CDS's phase is not 0 so we first need to fix that.
	# also the first/last is used later
	# these are already ordered

	# however don't really do this if you have a 5' partial, it will de-align to a potential splice site.
	# but sometimes we have a 5' partial with a UTR (which is terrible).

	my ($first_cds_exon_number,$last_cds_exon_number) = &get_cds_ends($isoform);
	my @exons = $isoform->get_exons();
	
	if ($exons[$first_cds_exon_number]->{CDS_exon_obj} && $exons[$first_cds_exon_number]->{CDS_exon_obj}->{'phase'} > 0){
		warn "FIXING: $error_name_text has the first CDS with a non-zero phase. Removing 5'UTR and adjusting...\n";
		#adjust exon and CDS, strip 5' UTR
		my ( $first_cds_end5,  $first_cds_end3 )  = $exons[$first_cds_exon_number]->{CDS_exon_obj}->get_coords();
		if ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{'strand'} eq '+'){
			$exons[$first_cds_exon_number]->set_coords( $first_cds_end5 + ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{'phase'}), $first_cds_end3 );
 		}else{
			$exons[$first_cds_exon_number]->set_coords( $first_cds_end5 - ($exons[$first_cds_exon_number]->{CDS_exon_obj}->{'phase'}) , $first_cds_end3 );
		}
		my @exon_coords = $exons[$first_cds_exon_number]->get_coords();
		$exons[$first_cds_exon_number]->{CDS_exon_obj}->set_coords( @exon_coords );
		$exons[$first_cds_exon_number]->{CDS_exon_obj}->{'phase'} = int(0);
		$isoform->{mRNA_exon_objs} = 0;
	   	$isoform->{mRNA_exon_objs} = \@exons;
		$isoform->refine_gene_object(); 
	}
	# now we need to check for potential joins that are illegal because the internal CDS has an exon with different co-ordinates.
	for (my $e=1; $e < scalar(@exons); $e++){
		next if $e <= $first_cds_exon_number;
		last if $e >= $last_cds_exon_number;
		my $exon = $exons[$e];
		my $cds = $exons[$e]->get_CDS_obj();
		next if !$cds;
		
		my @cds_coords = $cds->get_coords();
		my @exon_coords = $exon->get_coords();
		unless ($cds_coords[0] == $exon_coords[0] && $cds_coords[1] == $exon_coords[1]){
			#warn Dumper $isoform;
			die "FATAL:  $error_name_text ($first_cds_exon_number,$last_cds_exon_number) has an internal CDS ($e) with different co-ordinates than its corresponding exon. This is likely to be a fused prediction, please split them before continuing\n";
		}
	}
	if ($do_phasing){
		eval { $isoform->set_CDS_phases( $genome_seq_ref ); }; #if ever needed: unless $preferences{'norephase'};
	}

	$isoform->refine_gene_object();
	$isoform->create_all_sequence_types( $genome_seq_ref, %params );
	return $isoform;
}
