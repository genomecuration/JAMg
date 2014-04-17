#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) createRepeatLandscape.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Juan Caballero     jcaballero@systemsbiology.org
##  Description:
##      Based on earlier work by Juan Caballero
##
#******************************************************************************
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************
#
# ChangeLog
#
#     $Log: createRepeatLandscape.pl,v $
#     Revision 1.7  2013/11/06 19:01:26  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

createRepeatLandscape.pl - Create a Repeat Landscape graph

=head1 SYNOPSIS

  createRepeatLandscape.pl [-version] -div *.divsum [-t "graph title"]
                           [-j] [-g # ]

=head1 DESCRIPTION

  Create a Repeat Landscape graph using the divergence summary data
  generated with the calcDivergenceFromAlign.pl script.


=head1 EXAMPLES

  Older ( pre 4.0.4 ) RepeatMasker dataset:

     ./calcDivergenceFromAlign.pl -s example.divsum -a example_with_div.align 
                                  example.align.gz
 
     This creates an additional file "example_with_div.align" which contains
     the added Kimura divergence field after each alignment.

     ./createRepeatLandscape.pl -div example.divsum > 
                                /home/user/public_html/example.html 


  On newer RepeatMasker dataset that already contains the Kimura divergence
  line following each interspersed repeat alignment:

     ./calcDivergenceFromAlign.pl -s example.divsum example.align.gz

     ./createRepeatLandscape.pl -div example.divsum > 
                                /home/user/public_html/example.html 


The options are:

=over 4

=item -version

Displays the version of the program

=item -g #

Manually set the genome size.  

=item -j 

Output javascript only and not a fully constructed HTML page.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2013 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>
Juan Caballero <jcaballero@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Data::Dumper;
use Getopt::Long;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name: open-4-0-5 $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name: open-4-0-5 $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name: open-4-0-5 $';
my $CVSIdTag   =
    '$Id: createRepeatLandscape.pl,v 1.7 2013/11/06 19:01:26 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-div=s',
                    '-t=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

if ( !-s $options{'div'} ) {
  print "\n\nMissing '-div' option or the file is missing\n";
  usage();
}

my $kout    = $options{'div'};
my $gsize   = 'genome.size';
my $out     = undef;
my $help    = undef;
my $verbose = undef;

# Main variables
my %size;
my %data;
my %rep;
my %div;
my %sum;

# order is IMPORTANT for plotting

my $graphLabels = [
                    [ 'Unknown',         '#999999' ],
                    [ 'DNA/Academ',      '#FF0000' ],
                    [ 'DNA/CMC',         '#FF200B' ],
                    [ 'DNA/Crypton',     '#FF3115' ],
                    [ 'DNA/Ginger',      '#FF3D1E' ],
                    [ 'DNA/Harbinger',   '#FF4825' ],
                    [ 'DNA/hAT',         '#FF512D' ],
                    [ 'DNA/Kolobok',     '#FF5A34' ],
                    [ 'DNA/Maverick',    '#FF623B' ],
                    [ 'DNA',             '#FF6A42' ],
                    [ 'DNA/Merlin',      '#FF7149' ],
                    [ 'DNA/MULE',        '#FF7850' ],
                    [ 'DNA/P',           '#FF7F57' ],
                    [ 'DNA/PiggyBac',    '#FF865E' ],
                    [ 'DNA/Sola',        '#FF8D65' ],
                    [ 'DNA/TcMar',       '#FF936C' ],
                    [ 'DNA/Transib',     '#FF9972' ],
                    [ 'DNA/Zator',       '#FF9F79' ],
                    [ 'RC/Helitron',     '#FF00FF' ],
                    [ 'LTR/DIRS',        '#006400' ],
                    [ 'LTR/Ngaro',       '#197214' ],
                    [ 'LTR/Pao',         '#2A8024' ],
                    [ 'LTR/Copia',       '#3A8F33' ],
                    [ 'LTR/Gypsy',       '#489E42' ],
                    [ 'LTR/ERVL',        '#57AE51' ],
                    [ 'LTR',             '#65BD61' ],
                    [ 'LTR/ERV1',        '#73CD70' ],
                    [ 'LTR/ERV',         '#81DD80' ],
                    [ 'LTR/ERVK',        '#90ED90' ],
                    [ 'LINE/L1',         '#00008B' ],
                    [ 'LINE',            '#251792' ],
                    [ 'LINE/RTE',        '#38299A' ],
                    [ 'LINE/CR1',        '#483AA2' ],
                    [ 'LINE/Rex-Babar',  '#554BAA' ],
                    [ 'LINE/L2',         '#625CB1' ],
                    [ 'LINE/Proto2',     '#6E6DB9' ],
                    [ 'LINE/LOA',        '#797EC0' ],
                    [ 'LINE/R1',         '#848FC8' ],
                    [ 'LINE/Jockey-I',   '#8FA1CF' ],
                    [ 'LINE/Dong-R4',    '#99B3D7' ],
                    [ 'LINE/R2',         '#A3C5DE' ],
                    [ 'LINE/Penelope',   '#ACD8E5' ],
                    [ 'Other',           '#4D4D4D' ],
                    [ 'Other/composite', '#7F7F7F' ],
                    [ 'SINE',            '#9F1FF0' ],
                    [ 'SINE/5S',         '#A637F1' ],
                    [ 'SINE/7SL',        '#AD49F2' ],
                    [ 'SINE/Alu',        '#B358F3' ],
                    [ 'SINE/tRNA',       '#B966F4' ],
                    [ 'SINE/tRNA-Alu',   '#BF74F4' ],
                    [ 'SINE/tRNA-RTE',   '#C481F5' ],
                    [ 'SINE/RTE',        '#C98EF6' ],
                    [ 'SINE/Deu',        '#CE9BF7' ],
                    [ 'SINE/tRNA-V',     '#D3A7F7' ],
                    [ 'SINE/MIR',        '#D7B4F8' ],
                    [ 'SINE/Sauria',     '#DFCDF9' ],
                    [ 'SINE/tRNA-7SL',   '#E2D9F9' ],
                    [ 'SINE/tRNA-CR1',   '#E5E5F9' ]
];

my $sp = $kout;
$sp = "hg19";
my $gs  = 0;
my $rec = 0;
my @rep = ();
print STDERR "Parsing $kout\n";
open F, "$kout" or die "cannot read $kout\n";
my $inCoverageSection = 0;
while ( <F> ) {
  chomp;
  my @a = split( /\s+/, $_ );

  if ( /^Coverage for each repeat/ ) {
    $inCoverageSection = 1;
    next;
  }
  if ( /^Genome Size\s+=\s+(\d+)/ ) {
    $gs = $1;
  }
  if ( $inCoverageSection && /^Div\s+\S/ ) {
    shift @a;
    foreach my $r ( @a ) {
      $r =~ s/\?//g;

      # Specific name changes requested by Arian
      $r =~ s#DNA/Chompy#DNA#;
      $r =~ s#DNA/CMC-Chapaev#DNA/CMC#;
      $r =~ s#DNA/CMC-Chapaev-3#DNA/CMC#;
      $r =~ s#DNA/CMC-EnSpm#DNA/CMC#;
      $r =~ s#DNA/CMC-Transib#DNA/CMC#;
      $r =~ s#DNA/En-Spm#DNA/CMC#;
      $r =~ s#DNA/PIF-Harbinger#DNA/Harbinger#;
      $r =~ s#DNA/PIF-ISL2EU#DNA/Harbinger#;
      $r =~ s#DNA/Tourist#DNA/Harbinger#;
      $r =~ s#DNA/AcHobo#DNA/hAT#;
      $r =~ s#DNA/Charlie#DNA/hAT#;
      $r =~ s#DNA/Chompy1#DNA/hAT#;
      $r =~ s#DNA/MER1_type#DNA/hAT#;
      $r =~ s#DNA/Tip100#DNA/hAT#;
      $r =~ s#DNA/hAT-Ac#DNA/hAT#;
      $r =~ s#DNA/hAT-Blackjack#DNA/hAT#;
      $r =~ s#DNA/hAT-Charlie#DNA/hAT#;
      $r =~ s#DNA/hAT-Tag1#DNA/hAT#;
      $r =~ s#DNA/hAT-Tip100#DNA/hAT#;
      $r =~ s#DNA/hAT-hATw#DNA/hAT#;
      $r =~ s#DNA/hAT-hobo#DNA/hAT#;
      $r =~ s#DNA/hAT_Tol2#DNA/hAT#;
      $r =~ s#DNA/Kolobok-IS4EU#DNA/Kolobok#;
      $r =~ s#DNA/Kolobok-T2#DNA/Kolobok#;
      $r =~ s#DNA/T2#DNA/Kolobok#;
      $r =~ s#DNA/MULE-MuDR#DNA/MULE#;
      $r =~ s#DNA/MULE-NOF#DNA/MULE#;
      $r =~ s#DNA/MuDR#DNA/MULE#;
      $r =~ s#DNA/piggyBac#DNA/PiggyBac#;
      $r =~ s#DNA/MER2_type#DNA/TcMar#;
      $r =~ s#DNA/Mariner#DNA/TcMar#;
      $r =~ s#DNA/Pogo#DNA/TcMar#;
      $r =~ s#DNA/Stowaway#DNA/TcMar#;
      $r =~ s#DNA/Tc1#DNA/TcMar#;
      $r =~ s#DNA/Tc2#DNA/TcMar#;
      $r =~ s#DNA/Tc4#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Fot1#DNA/TcMar#;
      $r =~ s#DNA/TcMar-ISRm11#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Mariner#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Pogo#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Tc1#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Tc2#DNA/TcMar#;
      $r =~ s#DNA/TcMar-Tigger#DNA/TcMar#;
      $r =~ s#DNA/Tigger#DNA/TcMar#;
      $r =~ s#DNA/Helitron#RC/Helitron#;
      $r =~ s#LTR/DIRS1#LTR/DIRS#;
      $r =~ s#LTR/ERV-Foamy#LTR/ERVL#;
      $r =~ s#LTR/ERV-Lenti#LTR/ERV#;
      $r =~ s#LTR/ERVL-MaLR#LTR/ERVL#;
      $r =~ s#LTR/Gypsy-Troyka#LTR/Gypsy#;
      $r =~ s#LTR/MaLR#LTR/ERVL#;
      $r =~ s#LINE/CR1-Zenon#LINE/CR1#;
      $r =~ s#LINE/I#LINE/Jockey-I#;
      $r =~ s#LINE/Jockey#LINE/Jockey-I#;
      $r =~ s#LINE/L1-Tx1#LINE/L1#;
      $r =~ s#LINE/R2-Hero#LINE/R2#;
      $r =~ s#LINE/RTE-BovB#LINE/RTE#;
      $r =~ s#LINE/RTE-RTE#LINE/RTE#;
      $r =~ s#LINE/RTE-X#LINE/RTE#;
      $r =~ s#LINE/telomeric#LINE/Jockey-I#;
      $r =~ s#SINE/B2#SINE/tRNA#;
      $r =~ s#SINE/B4#SINE/tRNA-Alu#;
      $r =~ s#SINE/BovA#SINE/tRNA-RTE#;
      $r =~ s#SINE/C#SINE/tRNA#;
      $r =~ s#SINE/CORE#SINE/RTE#;
      $r =~ s#SINE/ID#SINE/tRNA#;
      $r =~ s#SINE/Lys#SINE/tRNA#;
      $r =~ s#SINE/MERMAID#SINE/tRNA-V#;
      $r =~ s#SINE/RTE-BovB#SINE/RTE#;
      $r =~ s#SINE/tRNA-Glu#SINE/tRNA#;
      $r =~ s#SINE/tRNA-Lys#SINE/tRNA#;
      $r =~ s#SINE/V#SINE/tRNA-V#;
      $r =~ s#Unknown/Y-chromosome#Unknown#;
      $rep{$r} = 1;
      push @rep, $r;
    }
    $rec = 1;
  }
  elsif ( $rec == 1 ) {
    my $div = shift @a;
    $div{$div} = 1;
    for ( my $i = 0 ; $i <= $#a ; $i++ ) {
      my $per = 100 * $a[ $i ] / $gs;
      my $rep = $rep[ $i ];
      $data{$sp}{$rep}{$div} += $per;
      $sum{$sp}{$rep}        += $per;
    }
  }
}
close F;

my $colors = "";
my @sprep  = ();
my @div    = sort { $a <=> $b } keys %div;
foreach my $label ( @{$graphLabels} ) {
  if ( defined $sum{$sp}{ $label->[ 0 ] } ) {
    push @sprep, $label->[ 0 ] if ( $sum{$sp}{ $label->[ 0 ] } > 0 );
    $colors .= "\"$label->[1]\", ";
  }
}
$colors =~ s/, $//;

#
# div  class1 class2..
#  #    frac
#
my %pdata = ();
foreach my $div ( @div ) {
  next if $div > 50;
  my @tmp = ();
  foreach my $rep ( @sprep ) {
    my $per = 0;
    $per = $data{$sp}{$rep}{$div} if ( defined $data{$sp}{$rep}{$div} );
    push @tmp, $per;
  }
  $pdata{$div} = join( ', ', @tmp );
}

# Default parameters
my $title = 'Repeat Landscape';
$title = $options{'t'} if ( $options{'t'} );
my @names = @sprep;

# Make it safe for the web
$title =~ s/_/ /g;

print <<_HEADER_
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
  <head>
    <title>Repeat Landscape</title>
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Divergence');
_HEADER_
    ;

foreach my $name ( @names ) {
  print "          data.addColumn('number', '$name');\n";
}

print "          data.addRows(\[\n";
foreach my $div ( sort { $a <=> $b } keys %pdata ) {
  my $dd = $pdata{$div};
  print "          \[\'$div\', $dd\],\n";
}

print <<_TAIL_
        ]);
        var options = {
          animation: {duration: 10},
          title: '$title',
          hAxis: {title: 'Kimura substitution level (CpG adjusted)', showTextEvery: 5},
          vAxis: {title: 'percent of genome'},
          isStacked: 1,
          colors: [$colors]
        };
        var chart = new google.visualization.ColumnChart(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
      </script>
    </head>
    <body bgcolor="#ffffff" text="#000000" link="#525d76">
    <font size="+3">Repeat Landscape</font>
   <hr noshade="noshade" size="1">
   <div id="chart_div" style="width: 1000px; height: 600px;"></div>
   <p>&#169; RepeatMasker.org</p>  
   </body>
</html>
_TAIL_
    ;

1;
