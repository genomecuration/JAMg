#!/usr/bin/perl -w 

=pod

=head1 USAGE

 Run this to prepare the BioPerl NCBI taxonomy

 -d => NCBI taxonomy flatfile directory (def. /databases/ncbi_taxonomy)
 -force => Force indexes to be rebuilt

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
use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
my $force = int(0);
my $flatfile_dir=$RealBin.'/../databases/ncbi_taxonomy/';
GetOptions(
 	'd|flatfile_dir:s' => \$flatfile_dir,
	'force' => \$force,
);

my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
   print "Using $flatfile_dir as NCBI TaxDB directory\n";
   $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $flatfile_dir.'/nodes.dmp',
     -namesfile => $flatfile_dir.'/names.dmp', -directory => $flatfile_dir, -force => $force) if !-s $flatfile_dir.'/nodes' || $force;
   print "Building protein binaries\n";
   new_dict (in => "$flatfile_dir/gi_taxid_prot.dmp", out => "$flatfile_dir/gi_taxid_prot.bin") if !-s "$flatfile_dir/gi_taxid_prot.bin" || $force;
   print "Building nucleotide binaries\n";
   new_dict (in => "$flatfile_dir/gi_taxid_nucl.dmp", out => "$flatfile_dir/gi_taxid_nucl.bin") if !-s "$flatfile_dir/gi_taxid_nucl.bin" || $force;
   print "Done!\n";
}else{
   die "Cannot find correct taxonomy at $flatfile_dir\n";
}

