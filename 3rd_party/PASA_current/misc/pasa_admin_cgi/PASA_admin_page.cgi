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


my $MYSQL_DB_DIR = &Pasa_conf::getParam("MYSQLDATA");

my $PASA_ADMIN_DB = &Pasa_conf::getParam("PASA_ADMIN_DB");
my $PASA_ADMIN_EMAIL = &Pasa_conf::getParam("PASA_ADMIN_EMAIL");
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_rw_user = &Pasa_conf::getParam("MYSQL_RW_USER");
my $mysql_rw_password = &Pasa_conf::getParam("MYSQL_RW_PASSWORD");

my $pasa_db_relocate_flag = (&Pasa_conf::getParam("PASA_DB_RELOCATE") =~ /true/i) ? 1:0;

my $cgi = new CGI;
print $cgi->header('text/html');

print $cgi->start_html(-title=>"Launch PASA-2");

print <<_HEAD;
<center>

<img src="show_png.cgi?image=images/PASA_logo.jpg" align=middle />
<hr>

</center>

_HEAD

    ;

my %params = $cgi->Vars();


## Launch job requests:
if ($params{runPASAalignAssembly}) {
    &process_alignment_assembly_submission();

} 
elsif ($params{runPASAannotationComparison}) {
    &process_run_annot_comparison_submission();

}
elsif ($params{runPASAaltSplice}) {
    &process_run_alt_splice_analysis();
} 


## Input forms:
elsif ($params{AlignAssemblyForm}) {
    &print_align_assembly_form();

} 
elsif ($params{AlignmentComparisonForm}) {
    &print_alignment_comparison_form();
}
elsif ($params{AltSpliceRequestForm}) {
    &print_alt_splice_request_form();
    
} else {
    &print_option_menu();
}


print $cgi->end_html();

exit(0);



####
sub print_option_menu {
    print <<_EOMENU_;

    <table align=center width=800>
    <tr><td>
    The PASA pipeline: Includes <b>PASA</b> (<b>P</b>rogram to <b>A</b>ssemble <b>S</b>pliced <b>A</b>lignments) as well as the pipeline to generate transcript alignments, compare alignment assemblies to existing gene model annotations, update gene structure annotations based on transcript alignments, and automatically model new genes based on full-length cDNA containing alignment assemblies. This system as well as its original application is described  in <a href="http://nar.oupjournals.org/cgi/content/full/31/19/5654">Nucleic Acids Res. 2003 Oct 1;31(19):5654-66</a>
	</td></tr>
	<tr><td><p align=left>All current PASA databases are listed <a href="PASA_admin_contents.cgi">here</a>.</td></tr>
	<tr><td><p align=left><a href="PASA_queue.cgi">Visit the PASA queue</a></td></tr>
<tr><td>
    <p><b>Please choose from the following options:</b></p>
	<ul>
	    <li><b>Launch the <a href="PASA_admin_page.cgi?AlignAssemblyForm=1">transcript alignment assembly pipeline:</a></b>  Given a set of genome and transcript sequences, the programs BLAT and sim4 are used to align all the transcripts to the genome.  High quality overlapping transcript alignments are then assemblied into maximal alignment assemblies, using the novel PASA alignment assembly algorithm.</li><br>
	    


        <li><b>Perform an <a href="PASA_admin_page.cgi?AlignmentComparisonForm=1">annotation comparison:</a></b> Compare the gene structures inferred by alignment assemblies to the current version of the genome annotation; identify new genes, new splicing isoforms for existing genes, and conflicts between existing gene structures and the transcript alignments.</li><br>
	    
        <li><b>Analyze <a href="PASA_admin_page.cgi?AltSpliceRequestForm=1">alternative splicing:</a></b> Pairwise comparisons are performed between PASA alignment assemblies found clustered to the same locus.  Specific splicing variations are identified including alternative acceptors and donors, retained introns, skipped exons, and alternative terminal exons. </li><br>
        

       <li><b>Automatically update gene structure annotations:</b>  New genes, isoforms of genes, and updates to existing gene structures can be performed automatically based on alignment assemblies. 

        <p>You can run the following script from the command line:
        <p><font size=-1>
        \$ANNOT_DEVEL/PASA2/scripts/cDNA_annotation_updater.dbi -P \"SYBTIGR,Sybase,\$sybase_username,\$sybase_password,\$annot_db\" -M \"\$pasa_db:$mysql_server:access:access\"
        <p>*Note the use of commas as separators for -P and colons as separators for -M.  Examine the additional options via -h
        
        <p>*Also note protein sequences in the database are NOT updated.  You will need to run <b>\$EGC_SCRIPTS/update_seqs.pl</b> to update all protein and CDS sequences stored in the database.
        

        </font>
        
        <li><b>Loading transcript alignments into the Euk Annot Database:</b>
        Just run the following, using your sybase database, username, and password:
        <p><font size=-1>
        \$ANNOT_DEVEL/PASA2/TIGR-only/annot_db_assembly_loader.dbi
        -S "\$db:SYBTIGR:\$username:\$password"
        -M "\$pasa_db,$mysql_server,access,access"
        -a 
        <p>examine the additional options via -h
        </font>


        </ul>

        <p>For assistance, please contact the current PASA administrator (<a href="mailto:$PASA_ADMIN_EMAIL">$PASA_ADMIN_EMAIL</a>) for this service.</b></li>	
		

		
</td></tr>
</table>

_EOMENU_

    ;
}



####
sub print_align_assembly_form {

print <<_EOHTML_;

<center>
<h1>PASA Job Submission Portal</h1>

<p align=left><a href="#dataRetrieval">Helpful hints</a> on generating input data sets are found at the bottom of this page.
<p align=left>* <i>required</i></p>
<form action="PASA_admin_page.cgi" method=post >

<table>
<tr bgcolor="#b686e0"><td>*Name of PASA MySQL database to be created:<br><i>(must end in <b>_pasa</b>)</i></td><td><input type=text name=MYSQLDB size=40 value=\"\*_pasa\"/></td></tr>
<tr bgcolor="#b686e0"><td>*Working Directory: Existing directory to post output and MySQL db files</td><td><input type=text name=WORKDIR size=80 /></td></tr>
<tr bgcolor="#b686e0"><td>*Your email address:</td><td><input type=text name=EMAIL size=40 /></td></tr>
</table>

<table>
<tr bgcolor="#b686e0"><td colspan=2 bgcolor="#006666" ><font color="#FFFFFF">Filenames as found in the Working Directory above:</font></td></tr>
<tr bgcolor="#86b6e0"><td width=50% >*Genome sequence fasta db<br><i><font color='#ff0000'><b>use asmbl_id for accession</b></font>, if applicable</i></td><td><input type=text name=GENOMIC_DB size=40 /></td></tr>
<tr bgcolor="#86b6e0"><td>*Transcript sequence fasta db<br>(<input type=checkbox name=VECTOR_POLYA_TRIM value=1 checked />trim vector and polyA)</td><td><input type=text name=TRANSCRIPT_DB size=40 /></td></tr>


<tr bgcolor="#86b6e0"><td>List of accessions for Full-length transcripts<br> (<i>highly recommended if known</i>) </td><td><input type=text name=FL_ACCS size=40 /></td></tr>

<tr><td colspan=2 bgcolor="#006666" ><font color= "#FFFFFF" > PASA configurable parameters for validating transcript alignments<br><i>Defaults are set</i></font></td></tr>
<tr bgcolor="#89e5e0"><td>Maximum intron length: </td><td><input type=text name=MAX_INTRON_LENGTH value=10000 /></td></tr>
<tr bgcolor="#89e5e0"><td>Minimum % transcript length aligned to genome:</td><td><input type=text name=MIN_PERCENT_ALIGNED value=90 ></td></tr>
<tr bgcolor="#89e5e0"><td>Minimum % identity for transcript alignment:</td><td><input type=text name=MIN_AVG_PER_ID value=95 ></td></tr>

<tr><td colspan=2 bgcolor="#006666" ><font color= "#FFFFFF" > Alignment Programs To Use</font></td></tr>
<tr  bgcolor="#79ffb7" ><td>Genome Mapping and Alignment Utility:</td>
   <td> 
    <input type=radio name="MAPPING_UTILITY" value="gmap" checked>gmap</input>
    <input type=radio name="MAPPING_UTILITY" value="blat" >blat</input>
 </td></tr>
<tr bgcolor="#79ffb7"><td>Retry with sim4 when above alignment fails validation:</td><td> <input type=checkbox name="SIM4_CHASER" /></td></tr>
<tr><td colspan=2 bgcolor="#006666" ><font color= "#FFFFFF" > Misc Opts</font></td></tr>
<tr bgcolor="#79ffb7"><td>Forcefully invalidate single-exon ESTs <br><i>*option INVALIDATE_SINGLE_EXON_ESTS</i></td>
     <td><input type=checkbox value="INVALIDATE_SINGLE_EXON_ESTS" /></td></tr>

<tr bgcolor="#79ffb7"><td>Use the Splice-Graph Assembler<br>
          <font size=-1><i>It is recommended that you use this if the product is for general annotation purposes.
                           If your goal is to examine alternative splicing, use the original PASA assembler and deselect this option.</i></font><br>
                           <i>*option USE_SPLICE_GRAPH_ASSEMBLER</i></td>
         <td><input type=checkbox name="USE_SPLICE_GRAPH_ASSEMBLER" /></td></tr>


<tr><td colspan=2 align=center><input type=submit value="Launch PASA Alignment &amp; Assembly" /></td></tr>

<tr><td colspan=2>
<p>After launching the PASA pipeline, you will be notified by email when it has been completed.  You will also be pointed to a URL from which you can navigate the results.</p>
</td></tr>
</table>
<input type=hidden name=runPASAalignAssembly value=1 />
</form>

</center>
<p>&nbsp;<p>
<h2><a name="dataRetrieval">Preparing PASA Input Data:</a></h2>
<ul>
 <li>Retrieving genome sequences from an Euk annotation database:<br>
 <font color='#ff0000'>echo select asmbl_id, sequence from assembly | runsql.dbi -D \$db -p ~/.pwdfile | runsql_to_fasta.dbi > genome.db </font><br>&nbsp;
 <li>Retrieving EST sequences from GenBank:<br>
     Visit Entrez Nucleotide database<br>
     Sample query: <font color='#ff0000'>oryza[Organism] AND "gbdiv est"[Properties]</font>
</ul>


_EOHTML_


    ;

}


####
sub process_alignment_assembly_submission {
    my ($MYSQLDB, 
        $WORKDIR, 
        $EMAIL, 
        $GENOMIC_DB, 
        $TRANSCRIPT_DB, 
        $FL_ACCS, 
        $MAX_INTRON_LENGTH,
        $MIN_PERCENT_ALIGNED,
        $MIN_AVG_PER_ID,
        $VECTOR_POLYA_TRIM,
        $MAPPING_UTILITY,
        $SIM4_CHASER) 
        = ($params{MYSQLDB},
           $params{WORKDIR},
           $params{EMAIL},
           $params{GENOMIC_DB},
           $params{TRANSCRIPT_DB},
           $params{FL_ACCS},
           int($params{MAX_INTRON_LENGTH}),
           int($params{MIN_PERCENT_ALIGNED}),
           int($params{MIN_AVG_PER_ID}),
           $params{VECTOR_POLYA_TRIM},
           $params{MAPPING_UTILITY},
           $params{SIM4_CHASER});
    my $INVALIDATE_SINGLE_EXON_ESTS = ($params{INVALIDATE_SINGLE_EXON_ESTS}) ? 1:0;
    my $USE_SPLICE_GRAPH_ASSEMBLER = ($params{USE_SPLICE_GRAPH_ASSEMBLER}) ? 1:0;
    
	print "<pre>" . Dumper (\%params) . "</pre>\n";
	#return;
	

    # set boolean flags
    $SIM4_CHASER = ($SIM4_CHASER) ? 1:0;
	

    unless ($VECTOR_POLYA_TRIM) {
        $VECTOR_POLYA_TRIM = 0;
    }
    
    my $errors = "";
    if (! ($MYSQLDB =~ /\w/ && $WORKDIR && $EMAIL && $GENOMIC_DB && $TRANSCRIPT_DB 
           && $MAX_INTRON_LENGTH > 0 && $MIN_PERCENT_ALIGNED >= 0 && $MIN_AVG_PER_ID >= 0 ) ) {
        $errors = "The form is incomplete. Please make sure you specify values for all fields.\n";
    } 
    else {
        
        unless ($MYSQLDB =~ /_pasa$/) {
            $errors .= "The pasa mysql database name must end in _pasa\n";
        }
        
        ## form completed. check values:
        unless (-d $WORKDIR) {
            $errors .= "directory [$WORKDIR] cannot be located.  Please check that this directory is available.\n";
        }
        
        unless (-w $WORKDIR) {
            $errors .= "I cannot write to directory [$WORKDIR].  Please relax your permissions to 0777";
        }
        
        ## reconstruct full paths to other files in the workdir:
        $GENOMIC_DB = "$WORKDIR/$GENOMIC_DB";
        $TRANSCRIPT_DB = "$WORKDIR/$TRANSCRIPT_DB";
                
        unless (-s $GENOMIC_DB) {
            $errors .= "$GENOMIC_DB genomic db cannot be found.\n";
        }
        
        unless (-s $TRANSCRIPT_DB) {
            $errors .= "$TRANSCRIPT_DB transcript db cannot be found.\n";
        }
        
        if ($FL_ACCS) {
            $FL_ACCS = "$WORKDIR/$FL_ACCS";
        
        	if (! -s $FL_ACCS) {
            	$errors .= "[$FL_ACCS] list of full-length accessions cannot be found.\n";
        	}	
        }    
		
        unless ($MAX_INTRON_LENGTH > 0) {
            $errors .= "value ($MAX_INTRON_LENGTH) specified for maximum intron length is invalid.\n";
        }
        unless ($MIN_PERCENT_ALIGNED > 0) {
            $errors .= "value ($MIN_PERCENT_ALIGNED) specified for minimum percent of the transcript to be aligned is invalid.\n";
        }
        unless ($MIN_AVG_PER_ID > 0) {
            $errors .= "value ($MIN_AVG_PER_ID) specified for minimum average percent identity for trancript alignment is invalid.\n";
        } 
    }
    if ($errors) {
        print "<PRE><font color='#FF0000'>Submission failed.  Errors:</font>\n$errors";
    } else {
        ## process submission:
        print "<pre>\n";
        eval {
            
            ## Build entry into PASA_admin
            my ($dbproc) = &connect_to_db($mysql_server,$PASA_ADMIN_DB,$mysql_rw_user, $mysql_rw_password);
            
            ## populate PASA_database_info
            my $query = "insert PASA_database_info (pasa_db_name, workdir, genome_db, transcript_db, fl_accs, timestamp, trim_vector_polyA)"
                . " values (?,?,?,?,?, now(), ?)";
            &RunMod($dbproc, $query, $MYSQLDB, $WORKDIR, $GENOMIC_DB, $TRANSCRIPT_DB, $FL_ACCS, $VECTOR_POLYA_TRIM);
            
            my $query = "select LAST_INSERT_ID()";
            my $db_id = &very_first_result_sql($dbproc, $query);
            
            ## populate audit_info
            my $query = "insert audit_info (pasa_db_name, job_type, email) values (?,?,?)";
            &RunMod($dbproc, $query, $MYSQLDB, "alignAssembly",  $EMAIL);
            
            my $query = "select LAST_INSERT_ID()";
            my $job_id = &very_first_result_sql($dbproc, $query);
            
            ## populate Alignment Assembly Parameters:
            my $query = "insert alignment_assembly_params (job_id, max_intron_length, min_percent_trans_align, min_percent_identity, mapping_utility, sim4_chaser_flag, invalidate_single_exon_ests, use_splice_graph_assembler) values (?,?,?,?, ?,?, ?, ?)";
            &RunMod($dbproc, $query, $job_id, $MAX_INTRON_LENGTH, $MIN_PERCENT_ALIGNED, $MIN_AVG_PER_ID, $MAPPING_UTILITY, $SIM4_CHASER, $INVALIDATE_SINGLE_EXON_ESTS, $USE_SPLICE_GRAPH_ASSEMBLER);
            
        };
        if ($@) {
            print "<pre><font color='#FF0000'>Sorry, submission failed for the following reasons:</font>\n$@\nPlease contact Brian Haas (bhaas\@tigr.org) for assistance.\n";
            
        } else {
            print "<h2>Submission has succeeded.  You will soon receive an email confirmation.</h2>\n";
        }
    }
    
    
}



####
sub print_alignment_comparison_form {

print <<_EOANNOTCOMPAREFORM_
    
<center>
<h1>PASA Annotation Comparison</h1>
<form action="PASA_admin_page.cgi" method=post >
<table border=1 width=800 >
<tr bgcolor="#b686e0"><td>Name of existing PASA MySQL database:</td><td><input type=text name=MYSQLDB size=40 /></td></tr>
<tr bgcolor="#b686e0">
    <td>Name of existing TIGR annotation database:</td><td><input type=text name=ANNOT_DB_NAME size=15 />
 	
	Load Gene Annotations:<input type=checkbox name=LOAD_LATEST_ANNOTS checked />
	
    </td>
</tr>
<tr bgcolor="#b686e0"><td>Your email address:</td><td><input type=text name=EMAIL size=40 /></td></tr>
<tr bgcolor="#b686e0"><td>Translation Table</td>
   <td><select name=GENETIC_CODE >
      <option value=universal >Universal</option>
      <option value=Euplotes >Euplotes</option>
      <option value=Tetrahymena >Tetrahymena</option>
      <option value=Candida >Candida</option>
      <option value=Acetabularia >Acetabularia</option>
   </select>
   </td></tr>


<!-- Mapping Requirements -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Requirements for Mapping Assemblies to Genes <i>(Defaults are set)</i></font></td></tr>


<tr bgcolor="#89e5e0"><td>Minimum % overlap between a gene and a full-length alignment assembly:<br><i>*option MIN_PERCENT_OVERLAP</i> </td>
<td><input type=text name=MIN_PERCENT_OVERLAP value=50 /></td></tr>

<!-- Requirements of genes based on FL-assemblies -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Requirements of genes based on FL-assemblies</font></td></tr>

<tr bgcolor="#89e5e0"><td>Minimum % of alignment assembly tentative cDNA sequence to be protein coding:<br><i>*option MIN_PERCENT_PROT_CODING</i></td>
<td><input type=text name=MIN_PERCENT_PROT_CODING value=40 ></td></tr>

<tr bgcolor="#89e5e0"><td>Minimum size of an ORF-encoded protein:<br><i>MIN_FL_ORF_SIZE<br><font size=-1 > (if set to zero (default), then only the MIN_PERCENT_PROT_CODING applies.  Otherwise, those failing the MIN_PERCENT_PROT_CODING test but encode an ORF of at least MIN_FL_ORF_SIZE will pass.) </font> </i></td>
<td><input type=text name=MIN_FL_ORF_SIZE value=0 ></td></tr>

<!-- Requirements of all genes -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Requirements for Valid Gene Models </font></td></tr>

<tr bgcolor="#89e5e0"><td>Maximum number of UTR exons:<br><i>*option MAX_UTR_EXONS</i></td>
<td><input type=text name=MAX_UTR_EXONS value=2 ></td></tr>

<!-- Comparisons among new genes and old genes, and new tentative isoforms  -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Comparisons among existing genes and tentative updates, and among tentative isoforms</font></td></tr>

<tr><td colspan=2 bgcolor="#0000EE" align=center ><font color= "#FFFFFF" > Length Criteria</font></td></tr>
<tr bgcolor="#89e5e0"><td>Minimum % of the original protein length for comparisons involving FL-assemblies:<br><i>*option MIN_PERCENT_LENGTH_FL_COMPARE</i></td>
<td><input type=text name=MIN_PERCENT_LENGTH_FL_COMPARE value=70 ></td></tr>


<tr bgcolor="#89e5e0"><td>Minimum % of the original protein length for comparisons involving non-FL assemblies:<br><i>*option MIN_PERCENT_LENGTH_NONFL_COMPARE</i></td>
<td><input type=text name=MIN_PERCENT_LENGTH_NONFL_COMPARE value=70 ></td></tr>

<tr><td colspan=2 bgcolor="#0000EE" align=center ><font color= "#FFFFFF" > FASTA Alignment Criteria</font></td></tr>

<tr bgcolor="#89e5e0"><td>Minimum % identity allowed for protein pairwise comparisons using Fasta:<br><i>*option MIN_PERID_PROT_COMPARE<br>used when comparing sibling isoforms and when comparing tentative protein updates to existing protein annotations</td>
<td><input type=text name=MIN_PERID_PROT_COMPARE value=70 ></td></tr>

<tr bgcolor="#89e5e0"><td>Minimum % of the shorter protein length which must align to the updated protein or isoform<br><i>*option MIN_PERCENT_ALIGN_LENGTH<br>used when comparing sibling isoforms and when comparing tentative protein updates to existing protein annotations</td>
<td><input type=text name=MIN_PERCENT_ALIGN_LENGTH value=70 ></td></tr>


<!-- Non-homology -based gene model comparison methods -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Non-homology based gene model comparison controls </font></td></tr>

<tr bgcolor="#89e5e0"><td>minimum % of overlap required among the genome span of the ORF of each overlapping gene to allow gene merging.<br><i>*option MIN_PERCENT_OVERLAP_GENE_REPLACE</i></td>
<td><input type=text name=MIN_PERCENT_OVERLAP_GENE_REPLACE value=80 ></td></tr>


<tr bgcolor="#89e5e0"><td>Ignore alignment results (% identity, % length) in protein comparisons.  Only consider the genome span of each ORF.<br><i>Requires the above minimum percentage of overlap among the existing ORF genome span.<br><i>*option STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE</i></td>
<td><input type=checkbox  name=STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE></td></tr>


<!-- Miscellaneous Options -->

<tr><td colspan=2 bgcolor="#006666" align=center ><font color= "#FFFFFF" > Miscellaneous Options </font></td></tr>

<tr bgcolor="#89e5e0"><td>Trust the FL-status of a cDNA (beware, may supposed FL are not really that).<br><i>It is better to not trust the FL-status when comparing and integrating into existing gene structures, because many are not actually FL.<br>*option TRUST_FL_STATUS</i></td>
     <td><input type=checkbox  name=TRUST_FL_STATUS></td></tr>


<tr><td colspan=2 align=center bgcolor="#006666"align=center ><input type=submit value="Launch Annotation Comparison" /></td></tr>

<tr><td colspan=2>
<p>After launching the PASA pipeline, you will be notified by email when it has been completed.  You will also be pointed to a URL from which you can navigate the results.</p>
</td></tr>
</table>
<input type=hidden name=runPASAannotationComparison value=1 />
</form>



_EOANNOTCOMPAREFORM_


    ;

}


####
sub process_run_annot_comparison_submission {
    my $MYSQLDB = $params{MYSQLDB}; 
    my $ANNOT_DB_NAME = $params{ANNOT_DB_NAME};
    my $MIN_PERCENT_OVERLAP =  $params{MIN_PERCENT_OVERLAP};
    my $MIN_PERCENT_PROT_CODING =  $params{MIN_PERCENT_PROT_CODING};
    my $MIN_PERCENT_LENGTH_FL_COMPARE = $params{MIN_PERCENT_LENGTH_FL_COMPARE};
    my $MIN_PERCENT_LENGTH_NONFL_COMPARE = $params{MIN_PERCENT_LENGTH_NONFL_COMPARE};
    my $MIN_PERID_PROT_COMPARE = $params{MIN_PERID_PROT_COMPARE};
    my $MIN_FL_ORF_SIZE =  $params{MIN_FL_ORF_SIZE};
    my $MAX_UTR_EXONS = $params{MAX_UTR_EXONS};
    my $MIN_PERCENT_ALIGN_LENGTH = $params{MIN_PERCENT_ALIGN_LENGTH};
    my $MIN_PERCENT_OVERLAP_GENE_REPLACE =  $params{MIN_PERCENT_OVERLAP_GENE_REPLACE};
    my $STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE = $params{STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE};
    my $TRUST_FL_STATUS =  $params{TRUST_FL_STATUS};
    my $EMAIL = $params{EMAIL};
    my $GENETIC_CODE = $params{GENETIC_CODE};
    my $LOAD_LATEST_ANNOTS = ($params{LOAD_LATEST_ANNOTS}) ? 1 : 0;
    
    ## set boolean flags
    $STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE = ($STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE) ? 1:0;
    $TRUST_FL_STATUS = ($TRUST_FL_STATUS) ? 1:0;
    
    
    #print "<pre>" . Dumper (\%params);
    
    my $errors = "";
    
    unless ($MYSQLDB && $ANNOT_DB_NAME && $MIN_PERCENT_OVERLAP && $MIN_PERCENT_PROT_CODING && 
            $MIN_PERCENT_LENGTH_FL_COMPARE && $MIN_PERCENT_LENGTH_NONFL_COMPARE && $MIN_PERID_PROT_COMPARE &&
            ($MAX_UTR_EXONS >= 0) && $MIN_PERCENT_ALIGN_LENGTH && $MIN_PERCENT_OVERLAP_GENE_REPLACE) {
        $errors = "Sorry, some required fields were not populated or contain invalid data.\n";
    }
    unless ($MYSQLDB) {
        $errors .= "PASA MySQL database name of ($MYSQLDB) is invalid.\n";
    }
    unless ($ANNOT_DB_NAME) {
        $errors .= "TIGR annotation database name ($ANNOT_DB_NAME) is invalid.\n";
    }
    unless ($MIN_PERCENT_OVERLAP) {
        $errors .= "Value ($MIN_PERCENT_OVERLAP) for minimum percent overlap is invalid.\n";
    }
    unless ($MIN_PERCENT_PROT_CODING) {
        $errors .= "Value ($MIN_PERCENT_PROT_CODING) for minimum percent protein coding is invalid.\n";
    }
    unless ($MIN_PERCENT_LENGTH_FL_COMPARE) {
        $errors .= "Value ($MIN_PERCENT_LENGTH_FL_COMPARE) for minimum percent length under FL-assembly comparisons is invalid.\n";
    }
    
    unless ($MIN_PERCENT_LENGTH_NONFL_COMPARE) {
        $errors .= "Value ($MIN_PERCENT_LENGTH_NONFL_COMPARE) for minimum percent length under non-FL-assembly comparisons is invalid.\n";
    }
    
    if ($MIN_FL_ORF_SIZE && $MIN_FL_ORF_SIZE < 0) {
        $errors .= "Value ($MIN_FL_ORF_SIZE) for minimum ORF length of FL-assemblies is invalid.\n";
    }
    
    unless ($MIN_PERID_PROT_COMPARE) {
        $errors .= "Value ($MIN_PERID_PROT_COMPARE) for minimum percent identity during protein sequence comparisons is invalid.\n";
    }
    
    unless ($MIN_PERCENT_ALIGN_LENGTH) {
        $errors .= "Value ($MIN_PERCENT_ALIGN_LENGTH) for minimum % align length is invalid.\n";
    }
    
    unless ($MIN_PERCENT_OVERLAP_GENE_REPLACE) {
        $errors .= "Value ($MIN_PERCENT_OVERLAP_GENE_REPLACE) for minimum % overlap of orf-genome span for gene merging and aggressive gene replacement is invalid.\n";
    }
    
    unless ($MAX_UTR_EXONS >= 0) {
        $errors . "Value ($MAX_UTR_EXONS) for max num utr exons is invalid. \n";
    }
    unless ($EMAIL =~ /\@/) {
        $errors .= "Value ($EMAIL) is an invalid email address. \n";
    }
    
    
    unless ($errors) {
        my $dbproc;
        # check to see if it is the same as any current entry.
        # if a current entry doesn't exist, insert it.
        eval {
            $dbproc = &connect_to_db($mysql_server,$PASA_ADMIN_DB,$mysql_rw_user, $mysql_rw_password);
            my $query = "select count(*) from PASA_database_info where pasa_db_name = ?";
            my $count = &very_first_result_sql($dbproc, $query, $MYSQLDB);
            if ($count) {
                ## insert it:
                my $query = "update PASA_database_info set annot_db_name = ? where pasa_db_name = ?";
                &RunMod($dbproc, $query, $ANNOT_DB_NAME, $MYSQLDB);
                
            } else {
                # no entry for MYSQL db:
                $errors .= "There doesn't appear to be a record for database ($MYSQLDB).  Are you sure a PASA alignment assembly has been performed using this database?\n";
            	die "no record of $MYSQLDB ";
			}
			
			# check to see if any genes have been loaded:
			unless ($LOAD_LATEST_ANNOTS) {
				my $dbproc2 = &connect_to_db($mysql_server,$MYSQLDB,$mysql_rw_user, $mysql_rw_password);
				my $query = "select count(*) from annotation_admin";
				my $annot_loaded_count = &very_first_result_sql($dbproc2, $query);
				unless ($annot_loaded_count) {
					$errors .= "There do not appear to be any existing versions of the gene structure annotations loaded into $MYSQLDB";
				  	die "<br>Go back and select to load the current gene annotations";
				}
			}
				
        };
        if ($@) {
            $errors .= "<br>Problem interacting with MySQL database: $@\n";
        } else {
            
            ## Enter job into queue:
            eval {
           
                ## post to audit_info
                my $query = "insert audit_info (pasa_db_name, job_type, email) values (?,?,?)";
                &RunMod($dbproc, $query, $MYSQLDB, "annotCompare", $EMAIL);
                
                my $query = "select LAST_INSERT_ID()";
                my $job_id = &very_first_result_sql($dbproc, $query);
                
                ## enter parameters in annot_compare_params table:
                $query = "insert annot_compare_params (job_id, min_percent_overlap, min_percent_prot_coding, min_perid_prot_compare, "
                    . "min_percent_length_fl_compare, min_percent_length_nonfl_compare, min_fl_orf_size, min_percent_align_length, "
                    . "min_percent_overlap_gene_replace, stomp_high_percentage_overlapping_gene, trust_fl_status, max_utr_exons, genetic_code, load_latest_annots) "
                    . "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?) ";
                &RunMod($dbproc, $query, $job_id, $MIN_PERCENT_OVERLAP, $MIN_PERCENT_PROT_CODING, $MIN_PERID_PROT_COMPARE,
                        $MIN_PERCENT_LENGTH_FL_COMPARE, $MIN_PERCENT_LENGTH_NONFL_COMPARE, $MIN_FL_ORF_SIZE, 
                        $MIN_PERCENT_ALIGN_LENGTH,
                        $MIN_PERCENT_OVERLAP_GENE_REPLACE, $STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE, 
                        $TRUST_FL_STATUS, $MAX_UTR_EXONS, $GENETIC_CODE, $LOAD_LATEST_ANNOTS);
                
            };
            if ($@) {
                $errors .= "Problem posting job to MySQL database $MYSQLDB: $@\n";
            }
            
        }
        if ($dbproc) {
            $dbproc->disconnect;
        }
    }
    
    
    if ($errors) {
        print "<pre><font color='#FF0000'>Submission of Annotation Comparison job failed.\nHere are the errors encountered:\n$errors\n</font>"; 
    } else {
        ## Successful submission.
        print "<h2><font color='#0000FF'>Submission of Annotation Comparison Request was Successful</h2>";
        print "You will soon receive an email confirmation.\n";
    }
}


####
sub print_alt_splice_request_form {
    
    print <<_EOF_ALT_SPLICE_FORM;

<center>    
<h1>Perform Analysis of Alternative Splicing</h1>
<form action="PASA_admin_page.cgi" method=post >
<table>
<tr><td>PASA database name:</td><td><input type=text size=40 name=MYSQLDB /></td></tr>
<tr><td>Your email address:</td><td><input type=text name=EMAIL size=40 /></td></tr>
<tr bgcolor="#006666"><td colspan=2 align=center> <input type=submit value="Launch Analysis of Alternative Splicing" /></td></tr>

</table>

<input type=hidden name=runPASAaltSplice value=1 />
</form>




_EOF_ALT_SPLICE_FORM



}



####
sub process_run_alt_splice_analysis {
    my $MYSQLDB = $params{MYSQLDB};
    my $EMAIL = $params{EMAIL};
    
    $MYSQLDB =~ s/\s//g;

    unless ($MYSQLDB =~ /\w/ && $EMAIL =~ /\w/) {
        print "Please go back and fill out the entire form.\n";
        return;
    }
    

    ## first, make sure this mysql database exists.

    eval {
        my $dbproc = &connect_to_db($mysql_server,$PASA_ADMIN_DB,$mysql_rw_user, $mysql_rw_password);
        my $query = "select count(*) from PASA_database_info where pasa_db_name = \"$MYSQLDB\"";
        my $count = &very_first_result_sql($dbproc, $query);
        unless ($count) {
            print "Sorry, there's no record of pasa database name [$MYSQLDB].  Did you mistype it?\n";
            return;
        }
        
        ## insert the job:
        
        eval { 
            
            my $query = "insert audit_info (pasa_db_name, job_type, email) values (?,?,?)";
            &RunMod($dbproc, $query, $MYSQLDB, "altSplice",  $EMAIL);
        };
        
        if ($@) {
            # rethrow exception
            die $@;
        }
        
        # request succeeded:
        print "Your request for an analysis of alternative splicing has succeeded.  You will receive an email once the job has been started.\n";
     
    };
    

    if ($@) {
        # error occurred.
        print "Sorry, your request could not be completed.  Please try again or contact PASA_ADMIN_EMAIL for assistance.\n";
        print "<font color='#FF0000'>Error: $@ </font>";
        return;
    }
    
    return;
}


