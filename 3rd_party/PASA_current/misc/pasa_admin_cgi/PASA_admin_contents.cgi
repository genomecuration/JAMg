#!/usr/local/bin/perl

use Pasa_init;
use Pasa_conf;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Mysql_connect;
use strict;
use DBI;
use Data::Dumper;

$|++;

my $db = &Pasa_conf::getParam("PASA_ADMIN_DB");
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my $base_url = &Pasa_conf::getParam("BASE_PASA_URL");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


my $cgi = new CGI();
print $cgi->header();
print $cgi->start_html(-title=>"PASA Admin contents");

my $query = "select pasa_db_name, annot_db_name, workdir, genome_db, transcript_db, fl_accs, timestamp, trim_vector_polyA from PASA_database_info order by timestamp";
my @results = &do_sql_2D ($dbproc, $query);

foreach my $result (@results) {
    my ($pasa_db_name, $annot_db_name, $workdir, $genome_db, $transcript_db, $fl_accs, $timestamp, $trim_vector_polyA) = @$result;
    
    print <<_EOTABLE;

    <table border=1>
	<tr><th>PASA database</th><td><a href="$base_url/status_report.cgi?db=$pasa_db_name" target=_new >$pasa_db_name</a></td></tr>
	<tr><th>TIGR annot db</th><td>$annot_db_name</td></tr>
	<tr><th>working directory</th><td>$workdir</td></tr>
	<tr><th>genome fasta file</th><td>$genome_db</td></tr>
	<tr><th>transcripts fasta file</th><td>$transcript_db</td></tr>
	<tr><th>FL accs</th><td>$fl_accs</td></tr>
	<tr><th>Date</th><td>$timestamp</td></tr>
	<tr><th>Vector/polyA trimming</th><td>$trim_vector_polyA</td></tr>
    </table>
    <p>

_EOTABLE


;
}



print $cgi->end_html();

$dbproc->disconnect;

exit;

			       
