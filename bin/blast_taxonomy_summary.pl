#!/usr/bin/perl -w 

=pod

=head1 USAGE

 -in or empty argument => input file
 -xml => Input is BLASTXML (otherwise BLAST is expected)
 -d => NCBI taxonomy flatfile directory (def. /databases/ncbi_taxonomy)
 -ncbi:s => is an NCBI blast (GI IDs) 'protein' or 'nucleotide'; otherwise Uniprot/Uniref is assumed

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

None known so far.

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
my $flatfile_dir=$RealBin.'/../databases/ncbi_taxonomy/';
my ($is_xml,$in,$verbose,$do_whole_taxonomy,$ncbi_dictionary,$is_ncbi);
GetOptions(
	'in=s'=>\$in,
	'top:i' => \$hit_limit,
	'ncbi:s' => \$is_ncbi,
	'all' => \$do_whole_taxonomy,
	'xml' => \$is_xml,
 	'd|flatfile_dir:s' => \$flatfile_dir,
	'verbose'	=> \$verbose,
);

$in=shift||warn ("Give a BLAST output file from UniProt\n") && pod2usage unless $in;
pod2usage " -ncbi must protein or nucleotide\n" unless (!$is_ncbi || $is_ncbi=~/prot/i || $is_ncbi=~/nuc/i);
my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
   warn "Using $flatfile_dir as NCBI TaxDB directory\n" if $verbose;
   $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile',-nodesfile => $flatfile_dir.'/nodes.dmp',-namesfile => $flatfile_dir.'/names.dmp',-directory => $flatfile_dir);
   if ($is_ncbi){
   	if ($is_ncbi=~/prot/i && !-s "$flatfile_dir/gi_taxid_prot.bin" && -s "$flatfile_dir/gi_taxid_prot.dmp"){
		print "Building protein binaries\n";
        	use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
		new_dict (in => "$flatfile_dir/gi_taxid_prot.dmp", out => "$flatfile_dir/gi_taxid_prot.bin");
	   }
	if ($is_ncbi=~/nuc/i && !-s "$flatfile_dir/gi_taxid_nucl.bin" && -s "$flatfile_dir/gi_taxid_nucl.dmp"){
		print "Building nucleotide binaries\n";
        	use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
		new_dict (in => "$flatfile_dir/gi_taxid_nucl.dmp", out => "$flatfile_dir/gi_taxid_nucl.bin");
	}
	$ncbi_dictionary = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$flatfile_dir/gi_taxid_prot.bin") if  $is_ncbi=~/prot/i;
        $ncbi_dictionary = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$flatfile_dir/gi_taxid_nucl.bin") if  $is_ncbi=~/nuc/i;
   }

}else{
   warn "Using remote Entrez query\n" if $verbose;
}

my $query_num=int(0);
if (!-s "$in.taxa_hub"){
print "Parsing $in";
print " as XML" if $is_xml;
print " using NCBI GI taxonomy" if $is_ncbi;
print "\n";
if ($is_xml){
	$query_num=`grep -c 'Iteration_query-ID' < $in`;chomp($query_num);
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
	$query_num=`grep -c '^Query=' < $in`;chomp($query_num);
	system("egrep -o 'Tax=.+[A-Za-z]\$|^Query=.+|^gi\|[0-9]+' $in > $in.taxa_hub");
}
}
my $number=`wc -l < $in.taxa_hub`||die ("File is empty?\n");chomp($number);
print "Building hash from $query_num sequences with $number taxa hits.\n";
print "I will now use the top $hit_limit hits only...\n";

my ( $counter,$flag,%hash);
my $notax=int(0);
open (IN, "$in.taxa_hub");
while (my $ln=<IN>){
        $counter++;
	if ($ln=~/^Query=/){
		my $hit=0;
		while (my $ln2=<IN>){
			if ($ln2=~/^Query=/){
				seek(IN, -length($ln2), 1);
				last;
			}
			if (!$is_ncbi && $ln2=~/^Tax=(\w+\s?\w*)/){
				my $tax=$1;
				chomp($tax);
				$tax=~s/\s+RepID.*//;
				$hash{$tax}++;
			}elsif($is_ncbi && $ln2=~/^gi=(\d+)/){
	        	        my $gi=$1;
				my $tax =  $ncbi_dictionary->get_taxid($gi);
				#die "gi $gi has no taxonomy. This is not possible...\n" if !$tax;
				if (!$tax){
					$notax++;
				}else{
			                $hash{$tax}++;
				}
			}
			$hit++;
			last if $hit == $hit_limit;
		}
	}
}
close IN;
print "\nFinished parsing...\nQuerying taxonomy...\n";
warn "Warning: Failed to get $notax taxonomy IDs from GIs. Number of successful retrievals was ".scalar(keys %hash).". If you think this is unreasonable your NCBI taxonomy database is probably out of date\n" if $notax;
open (OUT,">$in.taxa.tsv");
$counter=0;
foreach my $tax (sort {$a cmp $b} keys %hash){
        $counter++;
        if ( $counter =~ /0000$/ ) {
                print STDERR $counter."          \r";
        }
	my ($one,$two,$lineage,$genspe);
	$tax=~/^([1-9A-Z][0-9a-z]+)\s*([a-z]*)/;
	$one=$1;
	$two=$2;
	if ($two && $two!~/sp\.?\b/){
#print "Searching $tax for -g $one -s $two\n";
		($genspe,$lineage)= &get_species($one,$two);
	}elsif ($one){
#print "Searching $tax for -g $one\n";
		($genspe,$lineage) = &get_species($one);
	}
	chomp($lineage) if $lineage;
	$genspe=~s/\s+\d+\s*$//;
	$lineage="Lineage not found" if !$lineage;
	print OUT "$genspe\t$tax\t".$hash{$tax}."\t$lineage\n";
}
close OUT;
system ("sort -t '	' -nrk 3,3 -o $in.taxa.sort $in.taxa.tsv");
unlink("$in.taxa.tsv");
rename("$in.taxa.sort","$in.taxa.tsv");
system ("sed -i -e '1i#Species name\tNCBI TaxID\tContigs with this as top-$hit_limit hit\tLineage of Taxon' $in.taxa.tsv");


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
		         	warn "Could not find NCBI taxid for $one $two.\n";
		        	return 'Unknown';
			}
		}else{
	          	warn "Could not find NCBI taxid for $one.\n" if !$two;
	        	return 'Unknown';
		}
        }
my (%all_taxa,%all_taxa_sorted);
my $sorted_ref=\%all_taxa_sorted;
my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid) ;
if (!$taxon){
          	warn "Could not query NCBI for NCBI taxid $ncbi_taxid.\n";
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
	warn "Cannot find genus/species for NCBI tax ID $ncbi_taxid\n";
	return "Unknown";
}elsif(!$species){
	warn "Cannot find species for genus ($genus)\n";
	return "Unknown";
}elsif(!$genus){
	warn "Cannot find genus for species ($species)\n";
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
if (!$do_whole_taxonomy && $genus ne 'Unknown'){
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
                $lineage = "class:$class;order:$order;family:$family;genus:$genus;species:$species;common:$common;ncbi:$ncbi_taxid\n";
}elsif($genus eq 'Unknown' && !$do_whole_taxonomy){
        my $print="";
        foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
                $print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';' unless $all_taxa_sorted{$sort}{'name'} eq 'unknown';
        }
        $lineage = $print."common:$common;ncbi:$ncbi_taxid\n";


}else{
        my ($print);
        foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
                $print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';';
        }
        $lineage = $print."common:$common;ncbi:$ncbi_taxid\n";
}
	return ($genus.' '.$species,$lineage);
}

