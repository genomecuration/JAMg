#!/usr/bin/perl -w 

=pod

=head1 USAGE

 -in or empty argument => input file
 -top How many hits to use
 -genome_fasta

 -d => NCBI taxonomy flatfile directory (def. /databases/ncbi_taxonomy)

 -xml => Input is BLASTXML (otherwise BLAST is expected)
 -ncbi:s => is an NCBI blast (GI IDs) 'protein' or 'nucleotide'; otherwise Uniprot/Uniref is assumed

 -class

=head1 AUTHORS

 Alexie Papanicolaou 1 2

        1 Max Planck Institute for Chemical Ecology, Germany
        2 Centre for Ecology and Conservation, University of Exeter, UK
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 BioPerl Bio/DB/Taxonomy/flatfile.pm has an issue in the last few lines, where the DESTROY command
 deletes the indexes. You will need to comment out these unlink commands

=cut


use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");

my $hit_limit = 1;
my $tax_tag = 'OX';
my $flatfile_dir=$RealBin.'/../databases/ncbi_taxonomy/';
my ($is_xml,$in,$verbose,$do_whole_taxonomy,$ncbi_dictionary,$is_ncbi,$only_class,$genome_fasta);
GetOptions(
	'in=s'=>\$in,
	'top:i' => \$hit_limit,
	'ncbi:s' => \$is_ncbi,
	'all' => \$do_whole_taxonomy,
	'xml' => \$is_xml,
 	'd|flatfile_dir:s' => \$flatfile_dir,
	'verbose'	=> \$verbose,
	'class' => \$only_class,
	'tag_tax' => \$tax_tag,
	'genome_fasta:s' => \$genome_fasta
);


$in = shift||warn ("Give a BLAST output file from UniProt\n") && pod2usage unless $in;
pod2usage " -ncbi must protein or nucleotide\n" unless (!$is_ncbi || $is_ncbi=~/prot/i || $is_ncbi=~/nuc/i);

my $genome_hash = &process_fasta($genome_fasta)  if $genome_fasta;

my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
   system("$RealBin/prepare_ncbi_taxonomy.pl -d  $flatfile_dir") if !-s $flatfile_dir.'/nodes';
   warn "Using $flatfile_dir as NCBI TaxDB directory\n" if $verbose;
   $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $flatfile_dir, -nodesfile => $flatfile_dir.'/nodes.dmp', -namesfile=> $flatfile_dir.'/names.dmp',-force => 0);
   if ($is_ncbi){
	$ncbi_dictionary = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$flatfile_dir/gi_taxid_prot.bin") if  $is_ncbi=~/prot/i;
        $ncbi_dictionary = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$flatfile_dir/gi_taxid_nucl.bin") if  $is_ncbi=~/nuc/i;
   }

}else{
   warn "Using remote Entrez query\n" if $verbose;
}

if (!-s "$in.taxa_hub"){
print "Parsing $in";
print " as XML" if $is_xml;
print " using NCBI GI taxonomy" if $is_ncbi;
print "\n";
if ($is_xml){
	open (IN,$in) || die $!;
	open (OUT,">$in.taxa_hub");
	while (my $ln=<IN>){
		next unless $ln=~/^\s*<Hit_def/ || $ln=~/Iteration_query-def/ || $ln=~/<Hit_id>/;
		if (!$is_ncbi && $ln=~/(Tax=[A-Z][a-z]+\s?[a-z]+\s?[a-z]*)\s?[A-Z]?/){
			my $hit = $1;
			$hit=~s/ RepID.+//;
			print OUT $hit."\n";
		}elsif($ln=~/<Iteration_query-def>(\S+)/){
			print OUT "Query=".$1."\n";
		}elsif($is_ncbi && $ln=~/<Hit_id>gi\|(\d+)/){
			print OUT "gi=$1\n";
		}else{
			
		}
	}
	close IN;
	close OUT;
  }else{
	open (IN,$in) || die $!;
	open (OUT,">$in.taxa_hub");
	while (my $ln=<IN>){
		if ($ln=~/($tax_tag=[0-9]+)\s/){
			print OUT $1."\n";
		}elsif($ln=~/^(Query=.+)/){
			print OUT $1."\n";
		}elsif($ln=~/^(gi\|[0-9]+)\s/){
			print OUT $1."\n";
		}
		#system("grep -r -o '$tax_tag=[0-9]+\s|^Query=.+|^gi\\|[0-9]+' $in > $in.taxa_hub");
	}
	close IN;
	close OUT;
  }
}
print "I will now use the top $hit_limit hits only and the top hit for scaffolds\n";

my ( $counter,$flag,%ncbi_taxes,%ncbi_lineage);
my $notax=int(0);
my $orig_sep = $/;
$/ = 'Query= ';

open (OUT,">$in.taxa.tsv");
open (IN, "$in.taxa_hub");
my %scaffolds;

  while (my $record=<IN>){
        chomp($record); next unless $record;
	my @lines = split($orig_sep,$record);
	my $query = shift(@lines);
	$query=~s/\s+.+//;
	my $scaffold;
	if ($query =~/^(\S+)\.g\d+\.t\d+$/){
	   $scaffold = $1;
	   $scaffolds{$scaffold}{'proteins'}++;
	}else{
		die "Not an Augustus prediction...".$query;
	}
	my (%hash,$hit);
	foreach my $ln2 (@lines){
			if (!$is_ncbi && $ln2=~/^$tax_tag=(\d+)/){
				my $tax=$1;
				my $lineage = $ncbi_taxes{$tax} if $ncbi_taxes{$tax};
				if (!$lineage){
					$lineage = &get_species($tax) ;
					$ncbi_taxes{$tax} = &get_species($tax);
				}
				next unless $lineage;
				$ncbi_lineage{$lineage} = $tax if !$ncbi_lineage{$lineage};
				$hash{$lineage}++;
				if ($scaffold && !$hit){
					$scaffolds{$scaffold}{'lineage'}{$lineage}++;
				}
#				die Dumper ($tax, $lineage,\%hash,$query,$scaffold);	
			}
			$hit++;
			if ($hit == $hit_limit){
				print OUT $query."\t";
				foreach my $lineage (keys %hash){
					print OUT "\t".$hash{$lineage}."\t".$lineage."\t".$ncbi_lineage{$lineage};
				}
				print OUT "\n";
				last;
			}
	}
}
close IN;
close OUT;
#die Dumper \%ncbi_taxes;

$/ = $orig_sep;
if ($genome_fasta){
  open (OUT2,">$in.scaffolds.tsv2");
  print OUT2 "ID\tGAP_No\tScaffold_length\tProtein_No\tProteins_per_100kb\tResult_array\n";
  foreach my $scaffold (sort { $genome_hash->{$b}->{'length'} <=>$genome_hash->{$a}->{'length'}  } keys %{$genome_hash}){
	my $proteins = $scaffolds{$scaffold}{'proteins'} ? $scaffolds{$scaffold}{'proteins'} : int(0);
	
	my $print = $scaffold."\t" . $genome_hash->{$scaffold}->{'gaps'} ."\t" . $genome_hash->{$scaffold}->{'length'}."\t".$proteins
		   ."\t".sprintf( '%.2f', ($proteins / ($genome_hash->{$scaffold}->{'length'} / 100000) ) );
	foreach my $lineage (sort { $scaffolds{$scaffold}{'lineage'}{$b} <=> $scaffolds{$scaffold}{'lineage'}{$a}  } keys %{$scaffolds{$scaffold}{'lineage'}}){
		$print .= "\t".$scaffolds{$scaffold}{'lineage'}{$lineage}.':'.$lineage;
	}
	print OUT2 $print."\n";
 }
 close OUT2;
}

###########################
sub get_species(){
	my ($lineage,$ncbi_taxid);
	my $one = shift;
	return unless $one;
		my $two = shift;
	if ($one =~/^\d+$/){
		$ncbi_taxid = $one;
	}
	elsif ($two){
	    $ncbi_taxid = $taxondb->get_taxonid($one.' '.$two);
	}else{
	    $ncbi_taxid = $taxondb->get_taxonid($one);
	}
        if (!$ncbi_taxid){
		if ($two){
			$ncbi_taxid = $taxondb->get_taxonid($one);
			if (!$ncbi_taxid){
#		         	warn "Could not find NCBI taxid for $one $two.\n";
		        	return 'Unknown';
			}
		}else{
#	          	warn "Could not find NCBI taxid for $one.\n" if !$two;
	        	return 'Unknown';
		}
        }
my (%all_taxa,%all_taxa_sorted);
my $sorted_ref=\%all_taxa_sorted;
my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid) ;
if (!$taxon){
 #         	warn "Could not query NCBI for NCBI taxid $ncbi_taxid.\n";
        	return 'Unknown';
}
# we want class, order, family, genus, species, common name
my ($genus,$species,$class,$order,$family,$i);
 my $species_rank=$taxon->rank();
 my $species_name=$taxon->scientific_name();
 if ($species_rank eq 'genus'){
         $species = $species_name.' sp.';
         $genus=$species_name;
 }elsif($species_rank eq 'species'){
        $species = $species_name;
 }if (!$species){
        my $t=$taxon;
        while ($t=$t->ancestor()){
                $species_rank=$t->rank();
                next if !$species_rank;
                $species = $taxon->scientific_name() if ($species_rank eq 'species');
        }
        if (!$species){
                my $t=$taxon;
                while ($t=$t->ancestor()){
                        $species_rank=$t->rank();
                        next if !$species_rank;
                        $species = $taxon->scientific_name().' sp.' if ($species_rank eq 'genus');
                        $genus=$taxon->scientific_name();
                }
                if (!$species){
                        $species = 'Unknown';
                        $genus = 'Unknown';
                        $i++;$sorted_ref->{$i}={
                                'rank' => $species_rank,
                                'name' => $species_name,
                        };
                }
        }
 }
 $species =~s/^(\w+)\s//;
 $i++;$sorted_ref->{$i}={
        'rank' => 'species',
        'name' => $species,
 } if $species ne 'Unknown';
 $all_taxa{$species_rank}=$species;
 my $t=$taxon;
 unless ($genus && $genus eq 'Unknown'){
 while ($t=$t->ancestor()){
        my $rank=$t->rank();
        next if !$rank;
        my $name = $t->scientific_name();
        next if !$name;
        $i++;$sorted_ref->{$i}={
                'rank' => $rank,
                'name' => $name,
        };
        $all_taxa{$rank}=$name;
        if ($rank eq 'genus'){
                $genus=$name;
                last;
        }
 }
}
if (!$genus || !$species){
#	warn "Cannot find genus/species for NCBI tax ID $ncbi_taxid\n";
	return "Unknown";
}elsif(!$species){
#	warn "Cannot find species for genus ($genus)\n";
	return "Unknown";
}elsif(!$genus){
#	warn "Cannot find genus for species ($species)\n";
	return "Unknown";
}
 my $common="No common name";
 $common=$taxon->common_names() ? $taxon->common_names() : $genus.' '.$species if $genus ne 'Unknown';
 $t=$taxon;
 while ($t=$t->ancestor()){
        my $rank=$t->rank();
        next if !$rank;
        my $name = $t->scientific_name() ;
        next if !$name;
        $i++;
        $sorted_ref->{$i}={
                'rank' => $rank,
                'name' => $name,
        };
        $all_taxa{$rank}=$name;
        if ($rank eq 'class'){
                $class=$name;
                last if (!$do_whole_taxonomy && $genus ne 'Unknown');
        }elsif($rank eq 'order'){
                $order=$name;
        }elsif($rank eq 'family'){
                $family=$name;
        }
 }
if ($only_class){
	$lineage = $class;
}
elsif (!$do_whole_taxonomy && $genus ne 'Unknown'){
        if (!$class){
                # go up
                if ($all_taxa{'phylum'}){
                        $class=$all_taxa{'phylum'};
                }elsif ($all_taxa{'kingdom'}){
                        $class=$all_taxa{'kingdom'};
                }else{
                        $class='unknown';
                }
        }if (!$order){
                #go down one, then up one
                if ($all_taxa{'suborder'}){
                        $order=$all_taxa{'suborder'};
                }elsif ($all_taxa{'subclass'}){
                        $order=$all_taxa{'subclass'};
                }else{
                        $order='unknown';
                }
        }if (!$family){
                #go down one, then up one
                if ($all_taxa{'subfamily'}){
                        $family=$all_taxa{'subfamily'};
                }elsif ($all_taxa{'superfamily'}){
                        $family=$all_taxa{'superfamily'};
                }else{
                        $family='unknown';
                }
        }
                $lineage = "class:$class;order:$order;family:$family;genus:$genus;species:$species;common:$common;ncbi:$ncbi_taxid";
}elsif($genus eq 'Unknown' && !$do_whole_taxonomy){
        my $print="";
        foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
                $print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';' unless $all_taxa_sorted{$sort}{'name'} eq 'unknown';
        }
        $lineage = $print."common:$common;ncbi:$ncbi_taxid";
}else{
        my ($print);
        foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
                $print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';';
        }
        $lineage = $print."common:$common;ncbi:$ncbi_taxid";
}
	return $lineage;
}

####################
sub process_fasta($){
        my $file=shift;
        my $orig_sep = $/;
        $/=">";
        my %genome_sizes_hash;
        open (IN,$file) || die $!;
        while (my $rec=<IN>){
                chomp($rec);
                next unless $rec;
                my @lines=split("\n", $rec);
                my $id=shift(@lines);
                if ($id=~/^(\S+)\s/){$id=$1;}
                my $seq = join('',@lines);
                $seq=~s/\s+//g;
                next unless $id && $seq;
                $genome_sizes_hash{$id}{'length'} = length($seq);
                $genome_sizes_hash{$id}{'gaps'} = ( $seq =~ tr/N// );
        }
        close IN;
        $/=$orig_sep;
        return(\%genome_sizes_hash);
}


