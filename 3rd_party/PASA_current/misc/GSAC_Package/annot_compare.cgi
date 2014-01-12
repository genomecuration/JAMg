#!/util/bin/perl

use strict;
use warnings;
use CGI::Pretty ":standard";
use CGI::Carp qw(fatalsToBrowser);                   
use DB_File;
use lib ("PerlLib");
use Gene_cdna_image;
use String::CRC;
use Cwd;
use Storable qw (thaw);

$|++;

our $SEE = 0;
umask(0000);

my $javascript = "";
&init_globals(); # set javascript and anything else of interest.

unless ($ENV{WEBSERVER_TMP}) {
	$ENV{WEBSERVER_TMP} = "/tmp/images";

	if (! -d $ENV{WEBSERVER_TMP}) {
		mkdir ($ENV{WEBSERVER_TMP}) or die "Error, cannot create dir $ENV{WEBSERVER_TMP}";
	}

}

my $WRITE_IMAGE_ALWAYS = 1;

my $curr_dir = cwd;

## descriptions entered here (plain text or html) will show up in the hover boxes
my %descriptions = (
                    'COMPAT-end_agree' => '__TOKEN_A__ and __TOKEN_B__ structures are the same in their region of overlap and share either the very 5\' or 3\' terminal boundary',
                    'COMPLEX' => 'many-to-many gene mappings',
                    'ISOFORM_COMPAT-encaps' => '__TOKEN_A__ or __TOKEN_B__ model entirely consumes the other model structure and they are identical in the region of overlap',
                    'ISOFORM_COMPAT-staggered' => '__TOKEN_A__ and __TOKEN_B__ models have the same structure in their region of overlap, but have staggered boundaries',
                    'ISOFORM_DIFF' => '__TOKEN_A__ and __TOKEN_B__ gene models have different structures',
                    'ISOFORM_NOMAP' => 'Additional isoforms provided by __TOKEN_A__ or __TOKEN_B__ that lack counterparts.  Other isoforms of these genes must exist in either the ISOFORM_DIFF or ISOFORM_SAME categories',
                    'ISOFORM_SAME' => '__TOKEN_A__ and __TOKEN_B__ isoform structures encode identical CDS regions',
                    'MERGE' => 'A single __TOKEN_B__ gene maps to multiple __TOKEN_A__ genes',
                    'SPLIT' => 'Multiple __TOKEN_B__ genes map to a single __TOKEN_A__ gene',
                    '__TOKEN_A__-NOMAP' => '__TOKEN_A__-genes lacking a __TOKEN_B__ counterpart',
                    '__TOKEN_B__-NOMAP' => '__TOKEN_B__-genes lacking a __TOKEN_A__ counterpart',
                    'MAP_EXTREME_DIFF' => '__TOKEN_A__ and __TOKEN_B__ genes map positionally, but encode different reading frames',
    );


my $cgi = new CGI;


my $DEBUG = $cgi->param("DEBUG");
my $project = $cgi->param("project");
my $group = $cgi->param("group");
my $NUM_PER_PAGE = $cgi->param("NUM_PER_PAGE") || 10;

unless ($project && $group) {
    &generate_project_listing();
    exit(0);
}


my $comparison_file = "$group/$project/compare.out";
my $RAW = $cgi->param("raw");
unless (-s $comparison_file) {
    die "Error, cannot find $comparison_file";
}

if ($RAW) {
    &dump_raw_data();
    exit(0);
}

## parse the url templates
my $url_list_file = "$group/$project/URL_templates.txt";
my @URL_TEMPLATE_STRUCTS;
if (-s $url_list_file) {
    open (my $fh, $url_list_file) or die "Error, cannot open $url_list_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; } # comment line in template file
        unless (/\w/) { next; } # blank line not parsed
        my ($hyperlink_name, $url_template) = split (/\t/, $_, 2);
        push (@URL_TEMPLATE_STRUCTS, { 
                                        hyperlink => $hyperlink_name,
                                        url_template => $url_template,
                                    }
              );
    }
}


my $gene_obj_file = "$group/$project/gene_structures.inx";
unless (-s $gene_obj_file) {
    die "Error, cannot find gene structures file $gene_obj_file ";
}

my %index;
tie(%index, 'DB_File', $gene_obj_file, O_RDONLY, 0, $DB_BTREE);

unless (%index) {
    die "<p><font color='#FF0000'>Error resurrecting tied index</font>\n";
    
}


## opts that drive display
my $class = $cgi->param("class");
my $start = $cgi->param("start");
my $accession = $cgi->param("accession");

print $cgi->header;
if ($class && $start) {
    &generate_class_listing();
}
elsif ($accession) {
    &generate_accession_report();
}
else {
    &generate_entry_page();
}

exit(0);



#### 
sub generate_accession_report {
    my @captured_entries;

    print $cgi->start_html(-title => "Report for $accession",
                           -head  => style({type => 'text/css'},
                                           join('',<DATA>),
                                           ),
                           -script => $javascript,
                           );
    
    print "<h2>Report for $accession</h2>\n";
    
    open (my $fh, $comparison_file) or die "Error, cannot open $comparison_file";
    while (<$fh>) {
        chomp;
        my ($_class, $seqname, $model_listing) = split (/\t/);
        
        my $model_listing_copy = $model_listing;
        $model_listing_copy =~ s/:|;/ /g;
        
        if ($model_listing_copy =~ /\b$accession\b/) {
            push (@captured_entries, { seqname => $seqname,
                                       model_listing => $model_listing,
                                       class => $_class,
                                   } );
        }
        
    }
    close $fh;
    if (@captured_entries) {
        
        &generate_report_table(@captured_entries);
    }

    else {
        print "<font color='#FF0000'>Sorry, no entries found for accession:</font> [$accession]";
    }
    print $cgi->end_html();
    
}



####
sub generate_class_listing {
    
    print $cgi->start_html(
                           -title => "$project: Report for $class starting at $start",
                           -head  => style({type => 'text/css'},
                                           join('',<DATA>),
                                           ),
                           -script => $javascript,
                           
                           );
    
    print "<h2>$class</h2>\n";

    my $num_captured = 0;
    my $num_seen = 0;
    my @captured_entries;

    my $next_page_num = $start + $NUM_PER_PAGE;
    my $next_page_link = "<p><a href=\"annot_compare.cgi?group=$group&project=$project&class=$class&start=$next_page_num\">Next $NUM_PER_PAGE</a>\n"
        . "|"
            . "<a href=\"annot_compare.cgi?group=$group&project=$project\">_HOME_</a></p>\n";
    
    print $next_page_link;
    
    open (my $fh, $comparison_file) or die "Error, cannot open $comparison_file";
    while (<$fh>) {
        if (/\#/) { next;}
        chomp;
        my ($_class, $seqname, $model_listing) = split (/\t/);
        if ($_class eq $class) {
            $num_seen++;
            if ($num_seen >= $start) {
                
                push (@captured_entries, { seqname => $seqname,
                                           model_listing => $model_listing } );
                $num_captured++;
            }
        }
        
        if ($num_captured >= $NUM_PER_PAGE) {
            last;
        }
    }
    close $fh;
    
    &generate_report_table(@captured_entries);
    
    
    print $next_page_link;
    

    print $cgi->end_html();
}


####
sub generate_entry_page {
    
    print $cgi->start_html( 
                            -title => "Annotation Comparison to $project",
                            -head  => style({type => 'text/css'},
                                            join('',<DATA>),
                                            ),
                            -script => $javascript,
        );
    
    print "<h2>Annotation Comparison to: $project</h2>\n";

    my %class_counts;
    my %genes_classified;
    my %models_classified;

    my %types;
    
    open (my $fh, $comparison_file) or die "Error, cannot open $comparison_file";
    my $header = <$fh>;
    my ($gff_A_info, $token_A_info, $gff_B_info, $token_B_info) = split (/\t/, $header);
    my ($tokA_head, $token_A) = split (/\s+/, $token_A_info);
    my ($tokB_head, $token_B) = split (/\s+/, $token_B_info);
    
    
    while (<$fh>) {
        if (/^\#/) { next;}
        chomp;
        my ($class, $seqname, $model_listing) = split (/\t/);
        $class_counts{$class}++;
        foreach my $model_name (split (/,/, $model_listing)) {
            my ($gene, $model, $type) = split (/:+/, $model_name);
            $models_classified{$class}->{$type}->{$model} = 1;
            $genes_classified{$class}->{$type}->{$gene} = 1;
            $types{$type} = 1;
        }
    }
    close $fh;
    
    my @types = sort keys %types;
    print "<table border=1><tr>"
        . "<th>class</th><th>count</th>";
    foreach my $type (@types) {
        print "<th>$type</th>";
    }
    print "</tr>\n";
    
    my $row_toggle = 'odd';
    
    ## update description keys, replace tokens
    foreach my $class (keys %descriptions) {
        my $descr = $descriptions{$class};
        my $token_update_flag = 0;
        if ($class =~ /__TOKEN_A__/) {
            $class =~ s/__TOKEN_A__/$token_A/g;
            $token_update_flag = 1;
        }
        if ($class =~ /__TOKEN_B__/) {
            $class =~ s/__TOKEN_B__/$token_B/g;
            $token_update_flag = 1;
        }
        
        if ($token_update_flag) {
            $descriptions{$class} = $descr;
        }
        
    }
    


    foreach my $class (sort keys %class_counts) {
        print "<tr class='$row_toggle'>";
        $row_toggle = $row_toggle eq 'odd' ? 'even' : 'odd';
        
        ## check for a description
        if (! exists $descriptions{$class} ) {
            $descriptions{$class} = 'no description found';
        } else {
            ## perform token replacements
            $descriptions{$class} =~ s/__TOKEN_A__/$token_A/g;
            $descriptions{$class} =~ s/__TOKEN_B__/$token_B/g;
        }
        
        print "<td><a href=\"annot_compare.cgi?group=$group&project=$project&class=$class&start=1\">$class</a><div class='description'>$descriptions{$class}</div></td><td>$class_counts{$class}</td>";
        foreach my $type (@types) {
            my @genes = keys %{$genes_classified{$class}->{$type}};
            my $num_genes = scalar (@genes);
            
            my @models = keys %{$models_classified{$class}->{$type}};
            my $num_models = scalar (@models);
            
            print "<td>$num_genes genes, $num_models transcripts</td>";
        }
        print "</tr>\n";
    }
    print "</table>\n";
    
    &add_search_form();

    print "<p><a href=\"annot_compare.cgi?group=$group&project=$project&raw=1\">complete classification results in text format</a></p>\n";
    
    print $cgi->end_html();
}

####
sub add_search_form {
    
    print "<p>\n";
    print "<form method=get action=annot_compare.cgi>\n"
        . "<input type=hidden name=project value=$project>\n"
		. "<input type=hidden name=group value=$group >\n"
            . "<b>Identifier search</b>: <input type=text size=50 name=accession onChange=\"javascript:this.submit()\">\n"
                . "</form>\n";
    
}            



####
sub generate_report_table {
    my @entries = @_;
    
    unless ($start) {
        $start = 1;
    }
    
    my $count = $start;
    print "<table border=1>\n";
    

    foreach my $entry (@entries) {
        my $seqname = $entry->{seqname};
        my $model_listing = $entry->{model_listing};
        my $_class = $entry->{class};
        
        if ($_class) {
            print "<tr><td colspan=3><h2>$_class</h2></td></tr>\n";
        }
        
        my $crc = crc($model_listing, 32);
        
		my $image_filename = "$ENV{WEBSERVER_TMP}/$crc.png";
        
        my @gene_objs = &write_image($seqname, $model_listing, $image_filename);
        
        if ($SEE) {
            foreach my $gene_obj (@gene_objs) {
                print "<pre>" . $gene_obj->toString() . "</pre>\n";
            }
        }
                
        my $asmbl_id = $gene_objs[0]->{asmbl_id};
        
        my ($min_coord, $max_coord) = sort {$a<=>$b} $gene_objs[0]->get_coords();

        my %model_id_to_gene_obj;
        
        my %URL_tokens;

        foreach my $gene (@gene_objs) {
            
            my $model_id = $gene->{Model_feat_name};
            $model_id_to_gene_obj{$model_id}= $gene;
            my ($lend, $rend) = sort {$a<=>$b} $gene->get_coords();
            if ($lend < $min_coord) {
                $min_coord = $lend;
            }
            if ($rend > $max_coord) {
                $max_coord = $rend;
            }
        }
        if ($min_coord < 0) { $min_coord = 1;}
        
        my ($region_lend, $region_rend) = ($min_coord, $max_coord);
        $URL_tokens{"__REGION_LEND__"} = $region_lend;
        $URL_tokens{"__REGION_REND__"} = $region_rend;
        $URL_tokens{"__CONTIG_ID__"} = $seqname;

        $model_listing =~ s/,/, /g;
        
        
        ($min_coord, $max_coord) = sort {$a<=>$b} ($min_coord, $max_coord);
        
        my @models = split (/,/, $model_listing);
        
        print "<tr><td>$count</td>";
        
        my $div_count = 0;

        # the list of models and links:
        print "<td><ul class='model_links'>\n";
        foreach my $model (@models) {
            $div_count++;
            
            my ($tu_feat, $model_feat, $meta) = split (/:+/, $model);
        
            my $gene_obj = $model_id_to_gene_obj{$model_feat};
            
            my ($gene_end5, $gene_end3) = $gene_obj->get_coords();
            my $gene_id = $gene_obj->{TU_feat_name};
            
            $URL_tokens{"__GENE_ID__"} = $gene_id;
            $URL_tokens{"__MODEL_ID__"} = $model_feat;
            $URL_tokens{"__GENE_END5__"} = $gene_end5;
            $URL_tokens{"__GENE_END3__"} = $gene_end3;
            
            my @links = &get_links_with_replacements(\%URL_tokens);
            
            my $div_id = "menu.$count.$div_count";
            print "<li>" . &add_link_div($model, $div_id, \@links) . " </li>\n";
            
        }
        print "</ul></td>\n";
        
        ## add the image:
        print "<td><div class='shadowbox'><div class='shadow'><img src=\"show_png.cgi?image=$image_filename\" alt=\"$model_listing image\"></div></div></td>"
                                    . "</tr>\n";
        
        $count++;
    }
    
    print "</table>\n";

    &add_search_form();
    
        
}



####
sub write_image {
    
    my ($seqname, $model_listing, $image_filename) = @_;
    
    my @gene_objs;
    
    my @models = split (/,/, $model_listing);
    
    ## put working models first
    my @working = grep { /working/ } @models;
    my @other = grep { ! /working/ } @models;
    
    unless (@other) {
        #die "Error, no other in @models\n";
    }
    
    @models = (@working,@other);
    
    foreach my $model (@models) {
        my ($gene_id, $model_id, $ev_type) = split (/:+/, $model);
        my ($gene_obj) = &get_gene_obj($seqname, $model_id, $ev_type);
        if ($DEBUG) {
            print "<pre>" . $gene_obj->toString() . "</pre>\n";
        }
        
        push (@gene_objs, $gene_obj);
    }
    
    if ($WRITE_IMAGE_ALWAYS || ! -s $image_filename) {
        
        my $gene_image = Gene_cdna_image->new();
        
        my $image = $gene_image->create_image(@gene_objs);
    
		open (my $img_fh, ">$image_filename") or die $!;
        binmode $img_fh;
        print $img_fh $image->png();
        close $img_fh;
    }
    if ($DEBUG) {
        print "<pre>$image_filename\n</pre>\n";
    }
    return (@gene_objs);
    
}


####
sub get_gene_obj {
    my ($seqname, $feat_name, $ev_type) = @_;
    my $gene_obj_blob = $index{$feat_name} or die "Error, no gene_obj found for feat: $feat_name\n";
    my $gene_obj = thaw($gene_obj_blob);
    unless (ref $gene_obj) {
        die "Error, couldn't resurrect gene_obj for $feat_name, $ev_type";
    }
    
    $gene_obj->{Model_feat_name} = $feat_name;
    
    return ($gene_obj);
    
}



####
sub get_VectorBase_predictions {
    my ($seqname, @model_ids) = @_;
    
    my @gene_objs;

    my %model_data;

    my $transcript_id = "";
    my $gene_id = "";
    
    my %models_wanted;
    foreach my $model_id(@model_ids) {
        $models_wanted{$model_id}=1;
    }
    
    open (my $fh, "VectorBase_Aedes_aegypti_0.5.gff") or die $!;
    while (<$fh>) {
        chomp;
        my @x = split (/\t/);
        if ($x[0] ne $seqname) {
            next;
        }
        my ($seq_id, $feat_type, $end5, $end3, $orient, $gene_info) = ($x[0], $x[2], $x[3], $x[4], $x[6], $x[8]);
     
        
   
        if ($feat_type eq 'transcript') {
            $gene_info =~ /ID=(\S+); PARENT=(\S+)/;
            $transcript_id = $1;
            $gene_id = $2;
            $gene_id =~ s/;//;
        }
    
        unless ($models_wanted{$transcript_id}) {
            next;
        }
        
        if ($feat_type ne 'CDS') {
            next;
        }
        
        if ($orient eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
                       
        $model_data{$gene_id}->{$transcript_id}->{$end5} = $end3;
    }

    foreach my $gene_id (keys %model_data) {
        my @transcript_ids = keys %{$model_data{$gene_id}};
        
        my @curr_gene_objs;
        foreach my $transcript_id (@transcript_ids) {
       
            print "Adding: $gene_id, $transcript_id\n" if $SEE;
            my $coords_href = $model_data{$gene_id}->{$transcript_id};
            
            my $gene_obj = new Gene_obj();
            $gene_obj->populate_gene_obj($coords_href, $coords_href);
            
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{Model_feat_name} = $transcript_id;
            
            push (@curr_gene_objs, $gene_obj);
        }
        my $gene_obj = shift @curr_gene_objs;
        while (@curr_gene_objs) {
            my $isoform = shift @curr_gene_objs;
            $gene_obj->add_isoform($isoform);
        }
        $gene_obj->refine_gene_obj();
        push (@gene_objs, $gene_obj);
    }

    @gene_objs = sort {$a->{end5}<=>$b->{end5}} @gene_objs;
    
    return (@gene_objs);
}

####
sub get_broad_predictions {
    my ($seqname, @model_ids) = @_;
    
    # print "Getting broad preds: seqname($seqname), model_ids: (@model_ids)\n";
    
    my @gene_objs;

    my %model_data;

    my $transcript_id = "";
    my $gene_id = "";
    
    my %got_stop_codon;
    
    open (my $fh, "Aedes_V1_Manual_Transcripts.gtf") or die $!;
    while (<$fh>) {
        chomp;
        my @x = split (/\t/);
        if ($x[0] ne $seqname) {
            next;
        }
        my ($seq_id, $feat_type, $end5, $end3, $orient, $gene_info) = ($x[0], $x[2], $x[3], $x[4], $x[6], $x[8]);
             
        $gene_info =~ /gene_id \"(\S+)\"; transcript_id \"(\S+)\";/;
        $gene_id = $1 or die "no gene id";
        $transcript_id = $2 or die "no transcript id";
        
        if ($feat_type eq "stop_codon") {
            $got_stop_codon{$transcript_id}=1;
        }
        
        if ($feat_type ne 'CDS') {
            next;
        }
                        
        if ($orient eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
        
        $model_data{$gene_id}->{$transcript_id}->{$end5} = $end3;
    }
    
    my %models_wanted;
    foreach my $model_id (@model_ids) {
        $models_wanted{$model_id} = 1;
    }
    
    foreach my $gene_id (keys %model_data) {
        my @transcript_ids = keys %{$model_data{$gene_id}};
        
        my @curr_gene_objs;
        foreach my $transcript_id (@transcript_ids) {
    
            unless ($models_wanted{$transcript_id}) { next; }
            
            print "Adding: $gene_id, $transcript_id\n" if $SEE;
            my $coords_href = $model_data{$gene_id}->{$transcript_id};
            
            my $gene_obj = new Gene_obj();
            $gene_obj->populate_gene_obj($coords_href, $coords_href);
            
            $gene_obj->{TU_feat_name} = $gene_id . "::" . $transcript_id;
            $gene_obj->{Model_feat_name} = $gene_id . "::" . $transcript_id;
            
            if ($got_stop_codon{$transcript_id}) {
                
                ## add stop codon
                my $orient = $gene_obj->get_orientation();
                my @exons = reverse $gene_obj->get_exons();
                foreach my $exon (@exons) {
                    if (my $cds = $exon->get_CDS_obj()) {
                        if ($orient eq '+') {
                            $cds->{end3}+=3;
                        }
                        else {
                            $cds->{end3} -= 3;
                        }
                        last;
                    }
                }
                $gene_obj->refine_gene_object();
            }
            push (@curr_gene_objs, $gene_obj);
        }
        
        if (@curr_gene_objs) {
            my $gene_obj = shift @curr_gene_objs;
            while (@curr_gene_objs) {
                my $isoform = shift @curr_gene_objs;
                $gene_obj->add_isoform($isoform);
            }
            $gene_obj->refine_gene_obj();
            push (@gene_objs, $gene_obj);
        }
    }


    unless (@gene_objs) {
        die "Error, no gene objs found for @model_ids.\n";
    }
    @gene_objs = sort {$a->{end5}<=>$b->{end5}} @gene_objs;
    
    if ($SEE) {
        foreach my $gene_obj (@gene_objs) {
            print $gene_obj->toString();
        }
    } 
    
    return (@gene_objs);
}

####
sub dump_raw_data {
    print $cgi->header("text/plain");
    open (my $fh, $comparison_file) or die "Error, cannot open $comparison_file";
    while (<$fh>) {
        print;
    }
    close $fh;
}


#### 
sub generate_project_listing {

    print $cgi->header;
	
	my $project_list_file = "project_listing.txt";
    
	my $group = "";

	my %group_order;
	my $order = 0;

	my %group_to_eles;

	open (my $fh, $project_list_file) or die "Error, cannot find $project_list_file ";
    while (<$fh>) {
        if (/^\#/) { next; } # commented entry in file
        chomp;
        unless (/\w/) { next;}
        
		if (m|\[([^\]]+)\]|) {
			$group = $1;
		}
		else {
			my ($project, $description) = split (/\t/);
			
			unless ($group) { die "Error, group not set for $_ "; }

			unless (exists $group_order{$group}) {
				$order++;
				$group_order{$group} = $order;
			}
			
			push (@{$group_to_eles{$group}}, { 
				group => $group,
				project => $project,
				descr => $description,
			}
				  );
		}
	}
    close $fh;

    # build html
    
    print $cgi->start_html(-title => "Annotation Comparisons",
                           -head  => style({type => 'text/css'},
                                           join('',<DATA>),
                                           ),
                           -script => $javascript,
                           );
    
    print "<h1>Annotation Comparisons</h1>\n";
    
	
	foreach my $group (sort {$group_order{$a}<=>$group_order{$b}} keys %group_to_eles) {
		
		print "<div class='project_grouping'>\n";
		print "<h2 class='group_name'>$group</h2>\n";
		my @eles = @{$group_to_eles{$group}};
		
		print "<ul class='project_nav_list'>\n";
		
		foreach my $ele (@eles) {
			my ($group, $project, $descr) = ($ele->{group}, $ele->{project}, $ele->{descr});
			print "<li><a href=\"annot_compare.cgi?group=$group&project=$project\">$project</a>"
				. "<div class='description'>" . $descr . "</div>"
                . "</li>\n";
		}
		print "</ul>\n";
		
		print "</div>\n"; # end of project grouping.
	}
	
    print $cgi->end_html();

}


####
sub add_link_div {
    my ($model, $div_id, $links_aref ) = @_;
    

    my $link_text = "<a href=\"javascript:return(void);\" onClick=\"return clickreturnvalue()\" onMouseover=\"dropdownmenu(this, event, \'$div_id\')\">$model</a>\n";
    $link_text .= "<div id=\"$div_id\" class=\"anylinkcss\" style=\"width: 150px; background-color: lightyellow\">\n";

    foreach my $link (@$links_aref) {
        my $url = $link->{url};
        my $hyperlink = $link->{hyperlink};
        
        # gbrowse link:
        $link_text .= "<a href=\"$url\" target=_new >$hyperlink</a>\n";
    }
    
    $link_text .= "</div>\n";
    
    
    return ($link_text);
    
}

####
sub init_globals {
    $javascript = <<__EOJS__;


<!-- javascript from http://www.dynamicdrive.com/dynamicindex1/anylinkcss.htm -->

var disappeardelay=250  //menu disappear speed onMouseout (in miliseconds)
var enableanchorlink=0 //Enable or disable the anchor link when clicked on? (1=e, 0=d)
var hidemenu_onclick=1 //hide menu when user clicks within menu? (1=yes, 0=no)

/////No further editting needed

var ie5=document.all
var ns6=document.getElementById&&!document.all

function getposOffset(what, offsettype){
    var totaloffset=(offsettype=="left")? what.offsetLeft : what.offsetTop;
    var parentEl=what.offsetParent;
    while (parentEl!=null){
        totaloffset=(offsettype=="left")? totaloffset+parentEl.offsetLeft : totaloffset+parentEl.offsetTop;
        parentEl=parentEl.offsetParent;
    }
    return totaloffset;
}

function showhide(obj, e, visible, hidden){
    if (ie5||ns6)
        dropmenuobj.style.left=dropmenuobj.style.top=-500
    if (e.type=="click" && obj.visibility==hidden || e.type=="mouseover")
        obj.visibility=visible
    else if (e.type=="click")
        obj.visibility=hidden
}

function iecompattest(){
    return (document.compatMode && document.compatMode!="BackCompat")? document.documentElement : document.body
}

function clearbrowseredge(obj, whichedge){
    var edgeoffset=0
    if (whichedge=="rightedge"){
        var windowedge=ie5 && !window.opera? iecompattest().scrollLeft+iecompattest().clientWidth-15 : window.pageXOffset+window.innerWidth-15
        dropmenuobj.contentmeasure=dropmenuobj.offsetWidth
        if (windowedge-dropmenuobj.x < dropmenuobj.contentmeasure)
            edgeoffset=dropmenuobj.contentmeasure-obj.offsetWidth
    }
    else {
        var topedge=ie5 && !window.opera? iecompattest().scrollTop : window.pageYOffset
        var windowedge=ie5 && !window.opera? iecompattest().scrollTop+iecompattest().clientHeight-15 : window.pageYOffset+window.innerHeight-18
        dropmenuobj.contentmeasure=dropmenuobj.offsetHeight
        if (windowedge-dropmenuobj.y < dropmenuobj.contentmeasure) { //move up?
             edgeoffset=dropmenuobj.contentmeasure+obj.offsetHeight
             if ((dropmenuobj.y-topedge)<dropmenuobj.contentmeasure) //up no good either?
                 edgeoffset=dropmenuobj.y+obj.offsetHeight-topedge
        }
    }
       
    return edgeoffset
}

function dropdownmenu(obj, e, dropmenuID){
        
    if (window.event) event.cancelBubble=true
    else if (e.stopPropagation) e.stopPropagation()

    if (typeof dropmenuobj!="undefined") //hide previous menu
        dropmenuobj.style.visibility="hidden"
    clearhidemenu()
    if (ie5||ns6){
        obj.onmouseout=delayhidemenu
        dropmenuobj=document.getElementById(dropmenuID)
        if (hidemenu_onclick) dropmenuobj.onclick=function(){dropmenuobj.style.visibility='hidden'}
        dropmenuobj.onmouseover=clearhidemenu
        dropmenuobj.onmouseout=ie5? function(){ dynamichide(event)} : function(event){ dynamichide(event)}
        showhide(dropmenuobj.style, e, "visible", "hidden")
        /*    
        dropmenuobj.x=getposOffset(obj, "left")
        dropmenuobj.y=getposOffset(obj, "top")
        dropmenuobj.style.left=dropmenuobj.x-clearbrowseredge(obj, "rightedge")+"px"
        dropmenuobj.style.top=dropmenuobj.y-clearbrowseredge(obj, "bottomedge")+obj.offsetHeight+"px"
        */
    }
    return clickreturnvalue()
}

function clickreturnvalue(){
    if ((ie5||ns6) && !enableanchorlink) return false
    else return true
}

function contains_ns6(a, b) {
    while (b.parentNode)
        if ((b = b.parentNode) == a)
            return true;
    return false;
}

function dynamichide(e){
    if (ie5&&!dropmenuobj.contains(e.toElement))
        delayhidemenu()
    else if (ns6&&e.currentTarget!= e.relatedTarget&& !contains_ns6(e.currentTarget, e.relatedTarget))
        delayhidemenu()
}

function delayhidemenu(){
    delayhide=setTimeout("dropmenuobj.style.visibility='hidden'",disappeardelay)
}

function clearhidemenu(){
    if (typeof delayhide!="undefined")
        clearTimeout(delayhide)
}


__EOJS__

    return;

}

####
sub get_links_with_replacements {
    my ($URL_tokens_href) = @_;
    my @links;
    
    foreach my $url_template_struct (@URL_TEMPLATE_STRUCTS) {
        my ($hyperlink, $url_template) = ($url_template_struct->{hyperlink}, $url_template_struct->{url_template});
        
        foreach my $token (keys %$URL_tokens_href) {
            my $token_val = $URL_tokens_href->{$token};
            
            $url_template =~ s/$token/$token_val/g;
        }
        
        push (@links, {   
                          hyperlink => $hyperlink,
                          url => $url_template,
                      }
              );

    }

    return (@links);
}





################################################################################################
################################################################################################
### Fancy CSS to pretty up the page  ###########################################################
##################### by Joshua Orvis ##########################################################
################################################################################################
################################################################################################



__DATA__
/* this could be exported later if we want to keep it */


body {
    font-family: verdana,helvetica,sans-serif;
    font-size: 10pt;
}

a {
    color: rgb(25,25,255);
}

h1 {
    font-size: 120%;
    font-weight: bold;
    color: rgb(25,25,255);
    border-bottom: 1px dashed rgb(25,25,255);
    padding-bottom: 5px;
}

img {
    border: 1px solid rgb(165,165,165);
    padding: 3px;
    margin-bottom: 7px;
}

input, select {
    border: 1px solid rgb(0,0,0);
    padding-left: 3px;
}

table {
    border-style: none;
}

td, th {
    font-family: verdana,helvetica,sans-serif;
    font-size: 10pt;
    border-style: hidden;
    padding-left: 3px;
    padding-right: 3px;
}

td div.description, li div.description {
        display: none;
        position: absolute;
        border: 1px solid black;
        width: 35em;
        background-color: rgb(230,230,230);
        font-size: .9em;
        padding: .8em;
        margin-left: 7em;
        margin-top: 0.5em;
        opacity: 0.9;
}

td:hover div.description, li:hover div.description {
    display: block;
}

td div.description:hover, li.description:hover {
    display: none;
}

tr.odd td {
    background-color: rgb(225,225,225);
}

th {
    background-color: rgb(25,25,255);
    color: rgb(255,255,255);
    font-weight: bold;
}

.shadowbox {
	clear: both;
	float:left;
	background: url('/jorvis/aedes_annot_compare/shadow.gif') no-repeat bottom right;
	margin: 14px 0 0 17px !important;
	margin: 14px 0 0 8px;
}

.shadow {
    background: url('/jorvis/aedes_annot_compare/shadow2.png') no-repeat left top !important;
    background: url('/jorvis/aedes_annot_compare/shadow2.gif') no-repeat left top;
    float: left;
    padding: 0px 6px 6px 0px;
}

.shadowbox img, .shadowbox p {
    border: 1px solid #a9a9a9;
    margin: 0;
}


/* some additions by bhaas */

.project_nav_list li a { 
    text-decoration: none; 
}

.project_nav_list li, .model_links li, #other_aedes_links li {
    list-style: none;
    margin-top: 10px;
    padding: 0.25em;
    border-bottom: 1px solid gray;
  width: 180px;
}

.model_links li {
    width: 300px;
    font-size: 8pt;
    
}

.model_links a {
    text-decoration: none;
}


.project_grouping {
  background-color: #F3F3F3;
  padding: 0.25em;
  margin: 0.25em;
 }

.group_name  {
  position: relative;
  top: -10px;
  font-size: 100%;
  font-weight: bold;
  color: #660099;
}


/* from http://www.dynamicdrive.com/dynamicindex1/anylinkcss.htm */

.anylinkcss{
    position:absolute;
    left:200px;
    visibility: hidden;
    border:1px solid black;
    border-bottom-width: 0;
    font:normal 12px Verdana;
    line-height: 18px;
    z-index: 100;
    background-color: #E9FECB;
    width: 205px;
}

.anylinkcss a{
  width: 100%;
  display: block;
  text-indent: 3px;
  border-bottom: 1px solid black;
  padding: 1px 0;
  text-decoration: none;
  font-weight: bold;
  text-indent: 5px;
}

.anylinkcss a:hover{ /*hover background color*/
   background-color: black;
   color: white;
}
