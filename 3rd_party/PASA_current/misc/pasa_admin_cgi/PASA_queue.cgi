#!/usr/local/bin/perl

use Pasa_init;
use Pasa_conf;
use strict;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;

my $cgi = new CGI;
print $cgi->header('text/html');
my %params = $cgi->Vars();

my $MYSQL_DB_DIR = &Pasa_conf::getParam("MYSQLDATA");

my $PASA_ADMIN_DB = &Pasa_conf::getParam("PASA_ADMIN_DB");
my $PASA_ADMIN_EMAIL = &Pasa_conf::getParam("PASA_ADMIN_EMAIL");
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_rw_user = &Pasa_conf::getParam("MYSQL_RW_USER");
my $mysql_rw_password = &Pasa_conf::getParam("MYSQL_RW_PASSWORD");



my ($dbproc) = &connect_to_db($mysql_server,$PASA_ADMIN_DB,$mysql_rw_user, $mysql_rw_password);

if (my $job_id = $params{failed_job_id}) {
	&print_failed_job_error($job_id);
}
else {
	&report_queue();
}

print $cgi->end_html();

exit(0);

####
sub report_queue {

	my $query = "select job_id, pasa_db_name, is_started, start_time, is_finished, finish_time, is_success, job_type, email from audit_info order by job_id desc";
	my @results = &do_sql_2D($dbproc, $query);

	print "<table border=1 ><tr><th>PASA db</th><th>job type</th><th>start_time</th><th>finish_time</th><th>STATUS</th></tr>\n";

	foreach my $result (@results) {
		my ($job_id, $db_name, $is_started, $start_time, $is_finished, $finish_time, $is_success, $job_type, $email) = @$result;
	
		my $status;
		if ($start_time && $finish_time) {
			if ($is_success) {
				$status = "SUCCEEDED";
			}
			else {
				$status = "<a href=\"PASA_queue.cgi?failed_job_id=$job_id\">FAILED</a>";
			}
		}
		elsif ($start_time) {
			$status = "started";
		}
		else {
			$status = "pending";
		}
	
		print "<tr><th>$db_name</th>"
		. "<th>$job_type</th>"
	 	. "<th>$start_time</th>"
		. "<th>$finish_time</th>"
		. "<th>$status</th>"
		. "</tr>\n";
	}	

	print "</table>\n";


	return;
}



####
sub print_failed_job_error {
	my ($job_id) = @_;
	
	my $query = "select error_text from audit_info where job_id = $job_id";
	my $result = &very_first_result_sql($dbproc, $query);
	print "<pre>$result</pre>\n";
	
	return;
}
