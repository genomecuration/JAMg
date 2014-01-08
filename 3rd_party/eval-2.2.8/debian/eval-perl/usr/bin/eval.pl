#!/usr/bin/perl -w
use strict;
use lib '.';
use Tk;
use Tk::Pane;
use strict;
use Getopt::Std;
use Eval;
use GTF;

#constants
my $MIN_WIDTH = 440;
my $MIN_HEIGHT = 655;
my $GRAPHIC_MODE = 1;
my $TEXT_MODE = 2;
my $EVAL_FRAME_NUM = 0;
my $STATS_FRAME_NUM = 1;
my $FILTER_FRAME_NUM = 2;
my $GRAPH_FRAME_NUM = 3;
my $OVERLAP_FRAME_NUM = 4;
my $DIST_FRAME_NUM = 5;
my $ACTIVE_COLOR = 'black';
my $INACTIVE_COLOR = 'grey';
my $HOME = $ENV{HOME};
my @GENERAL_SKIP = qw(Correct
		      Specificity
		      Sensitivity
		      Matched);
my $MIN_DISP_LEN = 25;
my $MIN_VAL_LEN = 10;
my $FILTER_ALL = 1;
my $FILTER_ANY = 2;
my $USER = `whoami`;
chomp $USER;
my $GNUPLOT = $ENV{EVAL_GNUPLOT};
unless((defined($GNUPLOT)) && ($GNUPLOT ne "")){
    $GNUPLOT = "gnuplot";
}
my $GV = $ENV{EVAL_GV};
unless((defined($GV)) && ($GV ne "")){
    $GV = "gv";
}
my $HELP_PATH = $ENV{EVAL_HELP};
if((defined($HELP_PATH)) && ($HELP_PATH ne "")){
    if($HELP_PATH =~ /^(.+)\/$/){
	$HELP_PATH = $1;
    }
}
else{
    $0 =~ /^(.*)\/eval.pl$/;
    $HELP_PATH = "$1/help";
}


my @ALPHABET = qw(a b c d e f g h i j k l m n o p q r s t u v w x y z);

#global variables
my $X_Pos = 300;
my $Y_Pos = 100;
my $Options_File = "$HOME/.evalrc";
my @Main_Buttons;
my @Main_Frames;
my @Ann_Gtf_Lbs;
my @Pred_Gtf_Lbs;
my @Gtf_Objs;
my @Obj_Names;
my $Current_Frame;
my $Verbose = 0;
my $Really_Verbose = 0;
my $No_Seq = 0;
my $List_Mode = 0;
my $Mode = $GRAPHIC_MODE;
my $Precision = "%.2f";
my $Cwd = `pwd`;
my %General;
my %Display;
my %Graph_Resolution;
my $Next_Temp = 0;
my $Quick_Load = 0;
chomp $Cwd;

#parse command line options
use vars qw($opt_c $opt_g $opt_h $opt_l $opt_n $opt_q $opt_t $opt_v $opt_V);
getopts('cghlnqtvV');
my $usage = "usage: eval2.pl [-lnvVh] [gtf set 1] [gtf set 2] ...
Options:
   -c: Specify a the options file to load [default = ~/.evalrc] 
   -g: Load all files from the command line as single gtf files
   -l: Load all files from the command line as list files
   -n: No sequences.  Don't load conservation or fasta sequence.
   -q: Quick load.  Do not check gtf files for errors.
   -v: Verbose mode.  Outputs data to stdout.
   -V: Really verbose mode.  Standard verbose with any GTF format warnings included.
   -h: Help.  Displays this message.
";
$Verbose = ($opt_v || $opt_V);
$Really_Verbose = $opt_V;
$No_Seq = $opt_n;
if($opt_c){
    if(-e $opt_c){
	$Options_File = $opt_c;
    }
    else{
	error_func("File $opt_c does not exist.  Loading default file.");
    }
}
if($opt_l){
    $List_Mode = 1;
}
if($opt_g){
    $Mode = $GRAPHIC_MODE;
}
if($opt_t){
    $Mode = $TEXT_MODE;
}
if($opt_q){
    $Quick_Load = 1;
}
die $usage unless ((@ARGV >=2) or ($Mode == $GRAPHIC_MODE));
if($opt_h){
    print STDERR $usage;
    exit;
}

#main window
my $mw = MainWindow->new;
$mw->title("Eval");
$mw->geometry($MIN_WIDTH."x".$MIN_HEIGHT."+".$X_Pos."+".$Y_Pos);
$mw->minsize($MIN_WIDTH,$MIN_HEIGHT);
$mw->bind("<Destroy>",\&exit_func);
#Main_Frames
my $mw_menubar = $mw->Frame(-relief => 'ridge',
			    -borderwidth => 2)->pack(-side => 'top',
						     -fill => 'x');
my $mw_select_pane = $mw->Pane(-sticky => 'nswe')->pack(-side => 'top',
							-fill => 'x');
#make the above Pane scrolled to add a scrollbar to the button selection
#bar, incase it gets too long
my $mw_frame_select = $mw_select_pane->Frame(-relief => 'ridge',
					     -borderwidth => 2)->pack(-side => 'top',
								      -fill => 'both');
my $mw_main = $mw->Frame(-relief => 'ridge',
			 -borderwidth => 2)->pack(-side => 'top',
						  -expand => 1,
						  -fill => 'both');

$mw_main->packPropagate(0);

#create menus
$mw_menubar->Menubutton(-tearoff => 0,
			-text => "File",
			-underline => 0,
			-menuitems => [['command' => "Open",
					-underline => 0,
					-command  => \&open_func],
				       ['command' => "Save",
					-underline => 0,
					-command  => \&save_func],
				       ['command' => "Remove",
					-underline => 0,
					-command  => \&remove_func],
				       ['command' => "Exit",
					-underline => 1,
					-command  => \&exit_func]]
			)->pack(-side => 'left');

$mw_menubar->Menubutton(-tearoff => 0,
			-text => "Edit",
			-underline => 0,
			-menuitems => [['command' => "Options",
					-underline => 0,
					-command  => \&edit_options_func]]
			)->pack(-side => 'left');

$mw_menubar->Menubutton(-tearoff => 0,
			-text => "Help",
			-underline => 0,
			-menuitems => [['command' => "Help",
					-underline => 0,
					-command  => \&help_func],
				       ['command' => "About",
					-underline => 0,
					-command  => \&about_func]]
			)->pack(-side => 'right');



#create main frames
$Main_Frames[$EVAL_FRAME_NUM] = $mw_main->Frame()->pack(-side => 'top',
						   -expand => 1,
						   -fill => 'both');
$Main_Frames[$STATS_FRAME_NUM] = $mw_main->Frame();
$Main_Frames[$FILTER_FRAME_NUM] = $mw_main->Frame();
$Main_Frames[$GRAPH_FRAME_NUM] = $mw_main->Frame();
$Main_Frames[$OVERLAP_FRAME_NUM] = $mw_main->Frame();
$Main_Frames[$DIST_FRAME_NUM] = $mw_main->Frame();
$Main_Frames[$EVAL_FRAME_NUM]->Frame(-label => 'Evaluate Predictions',
				-borderwidth => 2,
				-relief => 'ridge')->pack(-side => 'top',
							  -fill => 'x');
$Main_Frames[$STATS_FRAME_NUM]->Frame(-label => 'Get Prediction Statistics',
				 -borderwidth => 2,
				 -relief => 'ridge')->pack(-side => 'top',
							   -fill => 'x');
$Main_Frames[$FILTER_FRAME_NUM]->Frame(-label => 'Filter Predictions',
				  -borderwidth => 2,
				  -relief => 'ridge')->pack(-side => 'top',
							  -fill => 'x');
$Main_Frames[$GRAPH_FRAME_NUM]->Frame(-label => 'Graph Predictions',
				 -borderwidth => 2,
				 -relief => 'ridge')->pack(-side => 'top',
							   -fill => 'x');
$Main_Frames[$OVERLAP_FRAME_NUM]->Frame(-label => 'Get Overlap Statistics',
					-borderwidth => 2,
					-relief => 'ridge')->pack(-side => 'top',
								  -fill => 'x');
$Main_Frames[$DIST_FRAME_NUM]->Frame(-label => 'Get Disribution',
					     -borderwidth => 2,
					     -relief => 'ridge')->pack(-side => 'top',
								       -fill => 'x');
$Current_Frame = $EVAL_FRAME_NUM;

#create main select buttons
$Main_Buttons[$EVAL_FRAME_NUM] = 
    $mw_frame_select->Button(-text => "Eval",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($EVAL_FRAME_NUM)});
$Main_Buttons[$STATS_FRAME_NUM] = 
    $mw_frame_select->Button(-text => "GenStats",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($STATS_FRAME_NUM)});
$Main_Buttons[$FILTER_FRAME_NUM] =
    $mw_frame_select->Button(-text => "Filter",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($FILTER_FRAME_NUM)});
$Main_Buttons[$GRAPH_FRAME_NUM] =
    $mw_frame_select->Button(-text => "Graph",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($GRAPH_FRAME_NUM)});
$Main_Buttons[$OVERLAP_FRAME_NUM] =
    $mw_frame_select->Button(-text => "Overlap",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($OVERLAP_FRAME_NUM)});
$Main_Buttons[$DIST_FRAME_NUM] =
    $mw_frame_select->Button(-text => "Dist",
			     -foreground => $INACTIVE_COLOR,
			     -relief => 'sunken',
			     -command => sub{switch_frames($DIST_FRAME_NUM)});

$Main_Buttons[$EVAL_FRAME_NUM]->grid($Main_Buttons[$STATS_FRAME_NUM],
				     $Main_Buttons[$FILTER_FRAME_NUM],
				     $Main_Buttons[$GRAPH_FRAME_NUM],
				     $Main_Buttons[$OVERLAP_FRAME_NUM],
				     $Main_Buttons[$DIST_FRAME_NUM],
				     -sticky => 'nsew', -ipadx => 5);
#initialize things here
init_hashes();
load_options_file();
initialize_frames();
switch_frames($EVAL_FRAME_NUM);

#load any files from command line
if($Verbose){
    print STDERR "Loading input files.\n";
}
if($Mode == $GRAPHIC_MODE){
    if(@ARGV > 0){
	foreach my $file (@ARGV){
	    load_func($file);
	}
    }
}

#Main Loop starts here
if($Verbose){
    print STDERR "Ready.\n";
}
MainLoop;
exit(0);

# init the General and Display hashes
sub init_hashes{
    my %order = Eval::get_list_struct();
    foreach my $level (@{$order{Levels}}){
	foreach my $type (@{$order{$level}{Type}}){
	    $Display{$level}{Type}{$type} = 1;
	    $General{$level}{Type}{$type} = 1;
	}
	foreach my $stat (@{$order{$level}{Stat}}){
	    $Display{$level}{Stat}{$stat} = 1;
	    $General{$level}{Stat}{$stat} = 1;
	    foreach my $skip (@GENERAL_SKIP){
		if($stat =~ /$skip/){
		    $General{$level}{Stat}{$stat} = 0;
		}
	    }
	}
    }
    foreach my $graph_type (Eval::get_graph_x_types()){
	my %default = get_default_graph_resolution();
	$Graph_Resolution{$graph_type} = \%default;
    }
}

# calls function to initialize all the frames 
sub initialize_frames{
    initialize_eval_frame();
    initialize_stats_frame();
    initialize_filter_frame();
    initialize_graph_frame();
    initialize_overlap_frame();
    initialize_distribution_frame();
}


# function called by the buttons on mw_frame_select to switch the frame currently
# being viewed
sub switch_frames{
    my ($new_frame) = @_;
    $Main_Frames[$Current_Frame]->packForget();
    $Main_Buttons[$Current_Frame]->configure(-relief => 'sunken',
					     -foreground => $INACTIVE_COLOR);
    $Main_Frames[$new_frame]->pack(-side => 'top',
			     -expand => 1,
			     -fill => 'both');
    $Main_Buttons[$new_frame]->configure(-relief => 'raised',
				   -foreground => $ACTIVE_COLOR);
    $Current_Frame = $new_frame;
}

#############################################################
#
#  Eval Frame Functions
# 
#############################################################
# initializes the eval preds frame
sub initialize_eval_frame{
    my ($frame,$alb,$plb) = make_ann_pred_frame($Main_Frames[$EVAL_FRAME_NUM]);
    $frame->configure(-borderwidth => 2,
		      -relief => 'ridge');
    $frame->pack(-side => 'top',
		    -fill => 'both',
		    -expand => 1);    
    my $button_frame = $Main_Frames[$EVAL_FRAME_NUM]->Frame(-borderwidth => 2,
							    -relief => 'ridge'
							    )->pack(-side => 'bottom',
								    -fill => 'x');
    $button_frame->Button(-text => "Run Eval",
			  -command => sub{eval_run_func($alb,$plb)}
			  )->pack(-side => 'right');
}

sub eval_run_func{
    my ($alb,$plb) = @_;    
    my @ann_index = $alb->curselection();
    my @pred_index = $plb->curselection();
    unless(($#ann_index >= 0) && ($#pred_index >= 0)){
	error_func("Must select at least one annotation and prediction.");
	return;
    }    
    my $ann = $Gtf_Objs[$ann_index[0]];
    my @preds;
    my @names = ($alb->get($ann_index[0]));
    foreach my $p_index (@pred_index){
	if($p_index == $ann_index[0]){
	    error_func("Cannot compare a gtf set to itself.");
	    $plb->selectionClear($p_index);
	    if($#pred_index == 0){
		return;
	    }
	}
	else{
	    push @preds, $Gtf_Objs[$p_index];
	    push @names, $Obj_Names[$p_index];
	}
    }
    my @data = Eval::evaluate($ann,\@preds,$Verbose);    
    #reduce number of digits after decimal point
    adjust_data_precision(\@data);
    display_eval_results(\@data,\@names);
}

sub display_eval_results{
    my ($data,$names) = @_;
    my $window = $mw->Toplevel(-title => "Evaluation Results");
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = $mw->height;
    my $width = $mw->height;
    my $report = "";
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $main = $window->Frame()->pack(-expand => 1,
				      -fill => 'both');
    my $main_pane = $main->Scrolled("Pane", 
				    -scrollbars => 'osoe',
				    -sticky => 'nswe')->pack(-fill => 'both',
							      -expand => 1);
    
    my $main_frame = $main_pane->Frame( 
					#-scrollbars => 'se',
				       -borderwidth => 2,
				       -relief => 'ridge')->pack(-side => 'top',
								 -expand => 1,
								 -fill => 'both');
    my $button_frame = $window->Frame(-borderwidth => 2,
				      -relief => 'ridge')->pack (-side =>'bottom',
								  -fill => 'x');
    $button_frame->Button(-text => 'Close',
			  -command => sub {$window->destroy}
			  )->pack(-side => 'right');
    $button_frame->Button(-text => 'Save',
			  -command => sub {save_eval_output($data,$names,$window)}
			  )->pack(-side => 'left');			  
    my $summary_frame = $main_frame->Frame(-label => 'Summary Statistics',
					   -borderwidth => 2,
					   -relief => 'ridge')->pack(-side =>'top',
								     -fill => 'both',
								     -expand => 1);
    my $general_frame = $main_frame->Frame(-label => 'General Statistics',
					   -borderwidth => 2,
					   -relief => 'ridge')->pack(-side =>'top',
								     -fill => 'both',
								     -expand => 1);
    my $details_frame = $main_frame->Frame(-label => 'Detailed Statistics',
					   -borderwidth => 2,
					   -relief => 'ridge')->pack(-side => 'top',
								     -fill => 'both',
								     -expand => 1);    
    #fill the summary stats frame
    for(my $i = 1; $i <= $#$data;$i++){
	my $frame = $summary_frame->Frame(-label => $$names[$i],
					  -borderwidth => 2,
					  -relief => 'ridge')->pack(-side => 'left',
								    -expand => 1,
								    -fill => 'both');
	my $label_frame = $frame->Frame()->pack(-side => 'left',
						-expand => 1,
						-fill => 'both');
	my $data_frame = $frame->Frame()->pack(-side => 'left',
					       -expand => 1,
					       -fill => 'both');
	$label_frame->Label(-text => "Gene Sensitivity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Gene}{All}{Consistent_Sensitivity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Gene Specificity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Gene}{All}{Consistent_Specificity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Transcript Sensitivity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Transcript}{All}{Consistent_Sensitivity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Transcript Specificity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Transcript}{All}{Consistent_Specificity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Exon Sensitivity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Exon}{All}{Correct_Sensitivity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Exon Specificity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Exon}{All}{Correct_Specificity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Nucleotide Sensitivity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Nuc}{All}{Correct_Sensitivity}
			   )->pack(-anchor => 'e');
	$label_frame->Label(-text => "Nucleotide Specificity")->pack(-anchor => 'w');
	$data_frame->Label(-text => $$data[$i]{Nuc}{All}{Correct_Specificity}
			   )->pack(-anchor => 'e');
    }
    fill_general_stats_frame($data,$names,$general_frame);
    #fill the detailed stats frame
    for(my $i = 1; $i <= $#$data;$i++){
	my $frame = $details_frame->Frame(-label => $$names[$i],
					  -borderwidth => 2,
					  -relief => 'ridge')->pack(-side => 'left',
								    -expand => 1,
								    -fill => 'both');
	my $label_frame = $frame->Frame()->pack(-side => 'left',
						-expand => 1,
						-fill => 'both');
	my $data_frame = $frame->Frame()->pack(-side => 'left',
					       -expand => 1,
					       -fill => 'both');	
	my %order = Eval::get_list_struct();	
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $label_frame->Label(-text => "$level_disp level"
				)->pack(-anchor => 'w');
	    $data_frame->Label(-text => "")->pack(-anchor => 'e');
	    foreach my $type (@{$order{$level}{Type}}){#keys %{$$data[$i]{$level}}){
		if($Display{$level}{Type}{$type} == 1){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    $label_frame->Label(-text => "   $type_disp"
					)->pack(-anchor => 'w');
		    $data_frame->Label(-text => "")->pack(-anchor => 'e');
		    foreach my $stat (@{$order{$level}{Stat}}){
			if($Display{$level}{Stat}{$stat} == 1){
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $label_frame->Label(-text => "      $stat_disp"
						)->pack(-anchor => 'w');
			    $data_frame->Label(-text => $$data[$i]{$level}{$type}{$stat}
					       )->pack(-anchor => 'e');
			}
		    }
		}
	    }
	}
    }
}

#fill the general stats frame
sub fill_general_stats_frame{
    my ($data,$names,$general_frame) = @_;
    for(my $i = 0; $i <= $#$data;$i++){
	my $frame = $general_frame->Frame(-label => $$names[$i],
					  -borderwidth => 2,
					  -relief => 'ridge')->pack(-side => 'left',
								    -expand => 1,
								    -fill => 'both');
	my $label_frame = $frame->Frame()->pack(-side => 'left',
						-expand => 1,
						-fill => 'both');
	my $data_frame = $frame->Frame()->pack(-side => 'left',
					       -expand => 1,
					       -fill => 'both');	
	my %order = Eval::get_general_list_struct();
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $label_frame->Label(-text => "$level_disp level"
				)->pack(-anchor => 'w');
	    $data_frame->Label(-text => "")->pack(-anchor => 'e');
	    foreach my $type (@{$order{$level}{Type}}){		
		if(($General{$level}{Type}{$type} == 1) &&
		   ($Display{$level}{Type}{$type} == 1)){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    $label_frame->Label(-text => "   $type_disp"
					)->pack(-anchor => 'w');
		    $data_frame->Label(-text => "")->pack(-anchor => 'e');
		    foreach my $stat (@{$order{$level}{Stat}}){
			if(($General{$level}{Stat}{$stat} == 1) &&
			   ($Display{$level}{Stat}{$stat} == 1) &&
			   ($stat !~ /Ann_/)){
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $label_frame->Label(-text => "      $stat_disp"
						)->pack(-anchor => 'w');
			    $data_frame->Label(-text => $$data[$i]{$level}{$type}{$stat}
					       )->pack(-anchor => 'e');
			}
		    }
		}
	    }
	}
    }
}

sub save_eval_output{
    my ($data,$names,$window) = @_;
    my $report = "";
    # add summary stats
    $report .= "\n**Summary Stats**\n";    
    $report .= "Annotation:\t$$names[0]\n";
    $report .= pad_string("Predictions:",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$names[$i]);
    }
    $report .= "\n";
    $report .= "\n".pad_string("Gene Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Gene}{All}{Consistent_Sensitivity});
    }
    $report .= "\n".pad_string("Gene Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Gene}{All}{Consistent_Specificity});
    }
    $report .= "\n".pad_string("Transcript Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Transcript}{All}{Consistent_Sensitivity});
    }
    $report .= "\n".pad_string("Transcript Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Transcript}{All}{Consistent_Specificity});
    }
    $report .= "\n".pad_string("Exon Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Exon}{All}{Correct_Sensitivity});
    }
    $report .= "\n".pad_string("Exon Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Exon}{All}{Correct_Specificity});
    }
    $report .= "\n".pad_string("Nucleotide Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Nuc}{All}{Correct_Sensitivity});
    }
    $report .= "\n".pad_string("Nucleotide Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Nuc}{All}{Correct_Specificity});
    }
    # add general stats
    $report .= get_general_stats_text($data,$names);
    # add detailed stats
    $report .= "\n\n**Detailed Stats**\n";
    $report .= "Annotation:\t$$names[0]\n";
    $report .= "Predictions:\t".pad_string("",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= pad_string($$names[$i]);
	if($i == $#$data){
	    $report .= "\n";
	}
	else{
	    $report .= "\t";
	}
    }
    if($#$data >= 0){
	my %order = Eval::get_list_struct();	
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $report .= "$level_disp\n";
	    foreach my $type (@{$order{$level}{Type}}){
		if($Display{$level}{Type}{$type} == 1){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    $report .= "\t$type_disp\n";
		    foreach my $stat (@{$order{$level}{Stat}}){
			if($Display{$level}{Stat}{$stat} == 1){			
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $report .= "\t\t".pad_string($stat_disp,$MIN_DISP_LEN);
			    for(my $i = 1; $i <= $#$data;$i++){
				$report .= "\t".
				    pad_string($$data[$i]{$level}{$type}{$stat});
			    }
			    $report .= "\n";
			}
		    }
		}
	    }
	}
    }
    my $filename = $window->getSaveFile(-defaultextension => ".eval",
					-filetypes        =>
					[['Eval Output','.eval'],
					 ['All Files','*']],
					-initialdir       => $Cwd,
					-initialfile      => "",
					-title            => "Save Eval Report As");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
	if(open(OUT,">$filename")){
	    print OUT $report;
	    message_func("File \"$filename\" saved.");
	    close(OUT);
	}
	else{
	    error_func("Could not open \"$filename\" for write.");
	}
    }
    else{
	error_func("No file selected.  File not saved");
    }
}

sub pad_string{
    my ($string,$min_len) = @_;
    unless(defined($min_len)){
	$min_len = $MIN_VAL_LEN;
    }
    my $len = length($string);
    for(my $i = $len;$i < $min_len;$i++){
	$string .= " ";
    }
    return $string;
}

sub get_general_stats_text{
    my ($data,$names) = @_;
    my $report = "";
    $report .= "\n\n**General Stats**\n";
    $report .= "Predictions:\n\t\t";
    for(my $i = 0;$i < $MIN_DISP_LEN;$i++){
	$report .= " ";
    }
    for(my $i = 0; $i <= $#$data;$i++){
	my $val = "$$names[$i]";
	my $len = length($val);
	for(my $i = $len;$i < $MIN_VAL_LEN;$i++){
	    $val .= " ";
	}
	$report .= "\t$val";
    }
    $report .= "\n";
    my %order = Eval::get_general_list_struct();
    if($#$data >= 0){
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $report .= "$level_disp\n";
	    foreach my $type (@{$order{$level}{Type}}){		
		if(($General{$level}{Type}{$type} == 1) &&
		   ($Display{$level}{Type}{$type} == 1)){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    my $len = length($type_disp);
		    for(my $i = $len;$i < $MIN_DISP_LEN;$i++){
			$type_disp .= " ";
		    }
		    $report .= "\t$type_disp\n";
		    foreach my $stat (@{$order{$level}{Stat}}){
			if(($General{$level}{Stat}{$stat} == 1) &&
			   ($Display{$level}{Stat}{$stat} == 1) &&
			   ($stat !~ /Ann_/)){
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $len = length($stat_disp);
			    for(my $i = $len;$i < $MIN_DISP_LEN;$i++){
				$stat_disp .= " ";
			    }
			    $report .= "\t\t$stat_disp";
			    for(my $i = 0; $i <= $#$data;$i++){
				my $val = $$data[$i]{$level}{$type}{$stat};
				$len = length($val);
				for(my $i = $len;$i < $MIN_VAL_LEN;$i++){
				    $val .= " ";
				}
				$report .= "\t$val";
			    }
			    $report .= "\n";
			}
		    }
		}
	    }
	}
    }
	return $report;
}

#############################################################
#
#  Stats Frame Functions
# 
#############################################################
# initializes the get pred stats frame
sub initialize_stats_frame{
    my ($frame,$plb) = make_pred_frame($Main_Frames[$STATS_FRAME_NUM]);
    $frame->configure(-borderwidth => 2,
		      -relief => 'ridge');
    $frame->pack(-side => 'top',
		 -fill => 'both',
		 -expand => 1);
    my $button_frame = $Main_Frames[$STATS_FRAME_NUM]->Frame(-borderwidth => 2,
						     -relief => 'ridge'
						     )->pack(-side => 'bottom',
							     -fill => 'x');
    $button_frame->Button(-text => "Get Stats",
			  -command => sub{get_stats_run_func($plb)}
			  )->pack(-side => 'right');
}


sub get_stats_run_func{
    my ($plb) = @_;
    my @pred_index = $plb->curselection();
    my @preds;
    my @names;
    foreach my $p_index (@pred_index){
	push @preds, $Gtf_Objs[$p_index];
	push @names, $Obj_Names[$p_index];
    }
    unless($#preds >= 0){
	error_func("Must select at least one prediction set.");
	return;
    }
    my @data = Eval::get_statistics(\@preds,$Verbose);
    #reduce number of digits after decimal point
    adjust_data_precision(\@data);
    display_stats_func(\@data,\@names);
}

sub display_stats_func{
    my ($data,$names) = @_;
    my $window = $mw->Toplevel(-title => "General Statistics");
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = $mw->height;
    my $width = $mw->width;
    my $report = "";
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $main = $window->Frame()->pack(-expand => 1,
				      -fill => 'both');
    
    my $main_pane = $main->Scrolled("Pane", 
				    -scrollbars => 'osoe',
				    -sticky => 'nswe')->pack(-fill => 'both',
							     -expand => 1);
    
    my $main_frame = $main_pane->Frame(-borderwidth => 2,
				       -relief => 'ridge')->pack(-side => 'top',
								 -expand => 1,
								 -fill => 'both');
    my $button_frame = $window->Frame(-borderwidth => 2,
				      -relief => 'ridge')->pack (-side =>'bottom',
								  -fill => 'x');
    $button_frame->Button(-text => 'Close',
			  -command => sub {$window->destroy}
			  )->pack(-side => 'right');
    $button_frame->Button(-text => 'Save',
			  -command => sub {save_stats_output($data,$names)}
			  )->pack(-side => 'left');			  
    fill_general_stats_frame($data,$names,$main_frame);
}

sub save_stats_output{
    my ($data,$names) = @_;
    my $report = get_general_stats_text($data,$names);
    my $filename = $mw->getSaveFile(-defaultextension => ".stats",
				    -filetypes        =>
				    [['Stats Output','.stats'],
				     ['All Files','*']],
				    -initialdir       => $Cwd,
				    -initialfile      => "",
				    -title            => 
				    "Save Statistics Report As");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
	if(open(OUT,">$filename")){
	    print OUT $report;
	    close(OUT);
	    message_func("File \"$filename\" saved.");
	}
	else{
	    error_func("Could not open \"$filename\" for write.");
	}
    }
    else{
	error_func("No file selected.  File not saved");
    }
    
}

#############################################################
#
#  Filter Frame Functions
# 
#############################################################
# initializes the filter preds frame
sub initialize_filter_frame{
    my $main_frame = $Main_Frames[$FILTER_FRAME_NUM]->Frame()->pack(-expand => 1,
								    -fill => 'both');
    my $ann_pred_frame = $main_frame->Frame()->pack(-side => 'top',
						    -fill => 'both',
						    -expand => 1);
    my $filter_frame = $main_frame->Frame();
    my ($frame,$alb,$plb) = make_ann_pred_frame($ann_pred_frame);
    $frame->configure(-borderwidth => 2,
		      -relief => 'ridge');
    $frame->pack(-side => 'top',
		 -fill => 'both',
		 -expand => 1);    
    my $ann_pred_button_frame 
	= $ann_pred_frame->Frame(-borderwidth => 2,
				 -relief => 'ridge'
				 )->pack(-side => 'bottom',
					 -fill => 'x');
    my $filter_filters_frame 
	= $filter_frame->Frame(-borderwidth => 2,
			       -relief => 'ridge')->pack(-side => 'top',
							 -fill => 'both',
							 -expand => 1);
    
    my $filter_button_frame 
	= $filter_frame->Frame(-borderwidth => 2,
			       -relief => 'ridge'
			       )->pack(-side => 'bottom',
				       -fill => 'x');
    my $filter_select_frame = $filter_filters_frame->Frame(-label => 
							   "Select A Filter To Add:",
							   -borderwidth => 0,
							   -relief => 'ridge',
							   )->pack(-padx => 20,
								   -pady => 10,
								   -side => 'top',
								   -fill => 'both',
								   -expand => 1);
    my $filter_key_frame = $filter_filters_frame->Frame(-label => "Filter Key:",
							-borderwidth => 0,
							-relief => 'ridge',
							)->pack(-padx => 20,
								-pady => 0,
								-side => 'top',
								-fill => 'x',
								-expand => 0);
    my $filter_string_frame = $filter_filters_frame->Frame(-label => 
							   "Enter Filter String:",
							   -borderwidth => 0,
							   -relief => 'ridge',
							   )->pack(-padx => 20,
								   -pady => 0,
								   -side => 'top',
								   -fill => 'x',
								   -expand => 0);
    my $filtered_name_frame = $filter_filters_frame->Frame(-label => 
							   "Enter string to append to ".
							   "new gtf name (optional):",
							   -borderwidth => 0,
							   -relief => 'ridge',
							   )->pack(-padx => 20,
								   -pady => 0,
								   -side => 'top',
								   -fill => 'x',
								   -expand => 0);
    my $filter_level_lb = $filter_select_frame->Scrolled("Listbox", 
							 -scrollbars => 'osoe',
							 -exportselection => 0,
							 -height => 6,
							 -selectmode => 'single'
							 )->pack(-side => 'left',
								 -fill => 'both',
								 -expand => 1,
								 -padx => 0,
								 -pady => 0);
    my %filter_types = Eval::get_filter_types();
    my @filter_type_lbs;
    my @type_map;
    foreach my $level (@{$filter_types{Levels}}){
	my $display = $level;
	$display =~ s/_/ /g;
	$filter_level_lb->insert('end',$display);
	my $lb = $filter_select_frame->Scrolled("Listbox", -scrollbars => 'osoe',
						-exportselection => 0,
						-height => 6,
						-selectmode => 'extended');
	foreach my $type (@{$filter_types{$level}}){
	    my $display = $type;
	    $display =~ s/_/ /g;
	    $lb->insert('end',$display);
	}
	$lb->selectionSet(0);
	push @filter_type_lbs, $lb;
	$type_map[$filter_level_lb->index('end')-1] = $#filter_type_lbs;
    }	
    $filter_level_lb->selectionSet(0);
    my $cur_type_lb = 0;
    $filter_type_lbs[$cur_type_lb]->pack(-side => 'right',
					 -fill => 'both',
					 -expand => 1,
					 -padx => 0,
					 -pady => 0);
    $filter_level_lb->bind("<Button>",sub{
	my ($lb) = @_;
	my @sel = $lb->curselection(); 
	if($cur_type_lb != $type_map[$sel[0]]){
	    my $a = 1;
	    $filter_type_lbs[$cur_type_lb]->packForget();
	    $cur_type_lb = $type_map[$sel[0]];
	    $filter_type_lbs[$cur_type_lb]->pack(-side => 'right',
						 -fill => 'both',
						 -expand => 1,
						 -padx => 0,
						 -pady => 0);
	}
    });
    my $filter_key_lb = $filter_key_frame->Scrolled("Listbox", -scrollbars => 'osoe',
						    -exportselection => 0,
						    -selectmode => 'extended'
						    )->pack(-side => 'top',
							    -fill => 'both',
							    -expand => 1);
    my $filter_text = "";
    my $filter_string_entry= $filter_string_frame->Entry(-textvariable => \$filter_text
							 )->pack(-side => 'bottom',
								 -fill => 'x');
    my $filtered_name = "";
    my $filtered_name_entry = $filtered_name_frame->Entry(-textvariable => 
							  \$filtered_name
							  )->pack(-side => 'bottom',
								  -fill => 'x');
    my %filter_keys;
    # Initialize Buttons
    $ann_pred_button_frame->Button(-text => "Select Filters ->",
				   -command => 
				   sub{
				       my @ann = $alb->curselection();
				       my @preds = $plb->curselection();
				       unless(($#ann >= 0) && ($#preds >= 0)){
					   error_func("Must select at least one ".
						      "annotation and prediction.");
					   return;
				       }
				       $ann_pred_frame->packForget();
				       $filter_frame->pack(-side => 'top',
							   -fill => 'both',
							   -expand => 1); 
				   }
				   )->pack(-side => 'right');
    $filter_button_frame->Button(-text => "<- Back",
				 -command => 
				 sub{
				     $filter_frame->packForget();
				     $ann_pred_frame->pack(-side => 'top',
							   -fill => 'both',
							   -expand => 1); 
				 }
				 )->pack(-side => 'left');
    $filter_button_frame->Button(-text => "Run",
				 -command => sub{
				     my @filter;
				     if(parse_filter_string(\%filter_keys,$filter_text,
							    \@filter)){
					 filter_run_func($alb,$plb,\@filter,
							 $filtered_name);
				     }
				 }
				 )->pack(-side => 'right');
    my %filter_in_use;
    $filter_button_frame->Button(-text => "Add Filter",
				 -command => sub{
                                     # add to filter_keys
				     my @indices = $filter_level_lb->curselection();
				     my $index = $indices[0];
				     my $level = $filter_level_lb->get($index);
				     my $filter_type_lb = $filter_type_lbs[$index];
				     @indices = $filter_type_lb->curselection();
				     foreach $index (@indices){
					 my $type = $filter_type_lb->get($index);
					 if(defined($filter_in_use{"$level$type"})){
					     next;
					 }
					 my $num = (keys %filter_keys);
					 my $letter = uc($ALPHABET[$num]);
					 $filter_key_lb->insert('end',"$letter - ".
								"$level $type");
					 $filter_key_lb->selectionSet($num);	
					 $filter_in_use{"$level$type"} = 1;
					 $level =~ s/ /_/g;
					 $type =~ s/ /_/g;
					 $filter_keys{$letter} = [$level,$type];
				     }
				 })->pack(-side => 'right');
    $filter_button_frame->Button(-text => "Remove Filter",
				 -command => sub{
				     $filter_text = "";
				     my @indices = $filter_key_lb->curselection();
				     my %map;
				     my @all = $filter_key_lb->get(0,'end');
				     my $next_check = 0;
				     my $next_index = 0;
				     $filter_key_lb->delete(0,'end');
				     for(my $i = 0;$i <= $#all;$i++){
					 my $text = $all[$i];
					 my @info = split /\s+/, $text;
					 my $level = $info[2];
					 my $type = join(" ",@info[3..$#info]);
					 if(($next_check <= $#indices) &&
					    ($indices[$next_check] == $i)){
					     delete($filter_in_use{"$level$type"});
					     $next_check++;
					 }
					 else{
					     my $letter = uc($ALPHABET[$next_index]);
					     $map{$letter} = [$info[2],$info[3]];
					     $filter_key_lb->insert('end',"$letter - ".
								    "$level $type");
					     $next_index++;
					 }
				     }
				     %filter_keys = %map;
				 })->pack(-side => 'right');
}

sub parse_filter_string{
    my ($filter_keys,$filter_text,$filter) = @_;
    $filter_text =~ s/ //g;
    $filter_text =~ s/AND/&/gi;
    $filter_text =~ s/OR/|/gi;
    while($filter_text =~ /&&/){
	$filter_text =~ s/&&/&/g;
    }
    while($filter_text =~ /\|\|/){
	$filter_text =~ s/\|\|/\|/g;
    }
    if($filter_text eq ''){
	error_func("Bad filter string format: Empty filter string.");
	return 0;
    }
    return parse_filter_helper($filter_keys,$filter_text,$filter);
}

sub parse_filter_helper{
    my ($keys,$text,$filter) = @_;
    my @chars = split //, $text;
    my $trim = 0;
    if($#chars == -1){
	error_func("Bad filter string format: Unmatched connector.");
	return 0;
    }
    elsif($chars[0] eq '('){
	my $in = 1;
	my $index = 1;
	while(($in > 0) && ($index <= $#chars)){
	    if($chars[$index] eq '('){
		$in++;
	    }
	    elsif($chars[$index] eq ')'){
		$in--;
	    }
	    $index++;
	}
	$index--;
	if($in == 0){
	    if($index == $#chars){
		return(parse_filter_helper($keys,join("",@chars[1..($index-1)]),
					   $filter));
	    }
	    else{
		my $start = $index;
		if($chars[$index+1] eq '&'){
		    $$filter[0] = 'And';
		    $start+=2;		
		}
		elsif($chars[$index+1] eq '|'){
		    $$filter[0] = 'Or';
		    $start+=2;		
		}
		else{
		    $$filter[0] = 'And';
		    $start = $index+1;
		}
		$$filter[1] = [];
		$$filter[2] = [];
		return(parse_filter_helper($keys,join("",@chars[0..$index]),
					   $$filter[1]) &&
		       parse_filter_helper($keys,join("",@chars[$start..$#chars]),
					   $$filter[2]));
	    }
	}
	else{
	    error_func("Bad filter string format: mismatched parentheses.");
	    return 0;
	}
    }
    else{
	my $f1_stop = -1;
	my $f2_start = -1;
	my $open = 0;
	my $close = 0;
	for(my $index = 1;$index < $#chars;$index++){
	    if($chars[$index] eq '('){
		$open++;
	    }
	    elsif($chars[$index] eq ')'){
		$close++;
	    }
	    if(($open - $close) == 0){
		if($chars[$index] eq '&'){
		    $$filter[0] = 'And';
		    $f1_stop = $index-1;
		    $f2_start = $index+1;
		    last;
		}
		elsif($chars[$index] eq '|'){
		    $$filter[0] = 'Or';
		    $f1_stop = $index-1;
		    $f2_start = $index+1;
		    last;
		}
	    }
	}
	# this string has a connector (ie AB&&CD, !(A&B)||C)
	if($f1_stop != -1){
	    $$filter[1] = [];
	    $$filter[2] = [];
	    return(parse_filter_helper($keys,join("",@chars[0..$f1_stop]),$$filter[1]) &&
		   parse_filter_helper($keys,join("",@chars[$f2_start..$#chars]),
				       $$filter[2]));
	}
	# this is a series of and connected filters (ie ABC, !ABC, !(ABC), !(AB)C)
	if($chars[0] eq '!'){
	    if($#chars == 0){
		error_func("Bad filter string format: unmatched negation.");
		return 0;
	    }
	    # this is a grouped negation (ie !(AB&C), !(AB)C)
	    if($chars[1] eq '('){
		my $in = 1;
		my $index = 2;
		while(($in > 0) && ($index <= $#chars)){
		    if($chars[$index] eq '('){
			$in++;
		    }
		    elsif($chars[$index] eq ')'){
			$in--;
		    }
		    $index++;
		}
		$index--;
		if($in == 0){
		    # this is a whole group (ie !(ABC))
		    if($index == $#chars){
			$$filter[0] = 'Not';
			$$filter[1] = [];
			return(parse_filter_helper($keys,join("",@chars[1..$#chars]),
						   $$filter[1]));
		    }
		    # this is a group negation and something more (ie !(AB)C, !(AB)!C)
		    else{
			$$filter[0] = 'And';
			$$filter[1] = [];
			$$filter[2] = [];
			return(parse_filter_helper($keys,join("",@chars[0..$index]),
						   $$filter[1]) &&
			       parse_filter_helper($keys,
						   join("",@chars[$index+1..$#chars]),
						   $$filter[2]));
		    }
		}
	    }
	    # this is a single filter negation
	    else{
		# this is asingle negation and more (ie !AB !A!B)
		if($#chars > 1){
		    $$filter[0] = 'And';
		    $$filter[1] = [];
		    $$filter[2] = [];
		    return(parse_filter_helper($keys,$chars[0..1],$$filter[1]) &&
			   parse_filter_helper($keys,join("",@chars[2..$#chars]),
					       $$filter[2]));
		}
		# this is a single negation (ie !A)
		else{
		    $$filter[0] = 'Not';
		    $$filter[1] = [];
		    return(parse_filter_helper($keys,join("",$chars[1]),
					       $$filter[1]));
		}
	    }
	}
	else{
	    # This is a string of letters (ie ABC)
	    if($#chars >= 1){
		$$filter[0] = 'And';
		$$filter[1] = [];
		$$filter[2] = [];
		return(parse_filter_helper($keys,$chars[0],$$filter[1]) &&
		       parse_filter_helper($keys,join("",@chars[1..$#chars]),
					   $$filter[2]));
	    }
	    # this is a single filter (ie A)
	    elsif($#chars == 0){
		unless(defined($$keys{$chars[0]})){
		    if(($chars[0] eq '&') ||
		       ($chars[0] eq '|') ||
		       ($chars[0] eq '!')){
			error_func("Bad filter string format: misplaced character, ".
				   "$chars[0].");
			return 0;
		    }		
		    else{
			error_func("Bad filter string format: bad character, ".
				   "$chars[0].");
		    }
		    return 0;
		}
		$$filter[0] = "Check";
		$$filter[1] = $$keys{$chars[0]}[0];
		$$filter[2] = $$keys{$chars[0]}[1];
		return 1;
	    }
	}
    }
    error_func("Bad filter string format.");
    return 0;
}

sub filter_run_func{
    my ($alb,$plb,$filter,$name) = @_;
    my @pred_index = $plb->curselection();
    my @preds;
    my @pred_names;
    foreach my $p_index (@pred_index){
	push @preds, $Gtf_Objs[$p_index];
	push @pred_names, $Obj_Names[$p_index];
    }
    my @ann_index = $alb->curselection();
    my $ann = $Gtf_Objs[$ann_index[0]];
    my $ann_name = $Obj_Names[$ann_index[0]];
    my @new_gtfs 
	= Eval::filter_predictions($ann,\@preds,$filter,$Verbose);
    if($name =~ /^\s+(\S.*)$/){
	$name = $1;
    }
    if($name =~ /^(.*\S)\s+$/){
	$name = $1;
    }
    if($name eq ""){
	$name = "filtered";
    }
    for(my $i = 0;$i <= $#new_gtfs;$i++){
	my $new_name = $pred_names[$i];
	if($new_name =~ /^(.*)\.gtf/){
	    $new_name = "$1.$name.gtf";
	}
	elsif($new_name =~ /^(.*)\.list/){
	    $new_name = "$1.$name.list";
	}
	else{
	    $new_name .= $name;
	}
	add_to_display($new_name,$new_gtfs[$i]);
    }
    message_func("Done Filtering");
}

#############################################################
#
#  Graph Frame Functions
# 
#############################################################
# initializes the graph preds frame
sub initialize_graph_frame{
    my $main_frame = $Main_Frames[$GRAPH_FRAME_NUM]->Frame()->pack(-expand => 1,
							      -fill => 'both');
    my $ann_pred_frame = $main_frame->Frame()->pack(-side => 'top',
						    -fill => 'both',
						    -expand => 1);
    my $graph_frame = $main_frame->Frame();
    my $select_frame = $main_frame->Frame();
    my ($frame,$alb,$plb) = make_ann_pred_frame($ann_pred_frame);
    $frame->configure(-borderwidth => 2,
		      -relief => 'ridge');
    $frame->pack(-side => 'top',
		 -fill => 'both',
		 -expand => 1);    
    my $ann_pred_button_frame 
	= $ann_pred_frame->Frame(-borderwidth => 2,
				 -relief => 'ridge'
				 )->pack(-side => 'bottom',
					 -fill => 'x');
    my $graph_graphs_frame 
	= $graph_frame->Frame(-borderwidth => 2,
			      -relief => 'ridge')->pack(-side => 'top',
							-fill => 'both',
							-expand => 1);
    my $graph_button_frame 
	= $graph_frame->Frame(-borderwidth => 2,
			      -relief => 'ridge'
			      )->pack(-side => 'bottom',
				      -fill => 'x');
    my $select_pred_frame = $select_frame->Frame(-borderwidth => 2,
						 -relief => 'ridge', 
						 -label => "Select Predictions:"
						 )->pack(-side => 'top',
							 -fill => 'both',
							 -expand => 1);

    my $select_graph_type_frame 
	= $select_frame->Frame(-borderwidth => 2,
			       -label => "Select Graphs to Display/View:",
			       -relief => 'ridge')->pack(-side => 'top',
							 -fill => 'both',
							 -expand => 1);
    my $select_graph_frame 
	= $select_frame->Frame(-borderwidth => 2,
			       -label => "Select Graphs to Display/View:",
			       -relief => 'ridge')->pack(-side => 'top',
							 -fill => 'both',
							 -expand => 1);
    my $select_button_frame 
	= $select_frame->Frame(-borderwidth => 2,
			       -relief => 'ridge'
			       )->pack(-side => 'bottom',
				       -fill => 'x');
    #add graph selection stuff here
    my $add_frm = $graph_graphs_frame->Frame(-label => "Select Graphs to Add:"
					     )->pack(-side => 'top',
						     -expand => 1,
						     -fill => 'x');
    my $y_frm = $add_frm->Frame()->pack(-side => 'top',-fill => 'x',
					-padx => 20);
    my $vs = $add_frm->Label(-text => "vs")->pack(-side => 'top');
    my $x_frm = $add_frm->Frame()->pack(-side => 'top',-fill => 'x');
    my $xlb = $x_frm->Scrolled("Listbox", -scrollbars => 'osoe',
			       -exportselection => 0,
			       -height => 5,
			       -selectmode => 'extended'
			       )->pack(-side => 'top',
				       -fill => 'both',
				       -expand => 1,
				       -anchor => 'w',
				       -padx => 20,
				       -pady => 0);
    my $x_level_lb = $y_frm->Scrolled("Listbox", -scrollbars => 'osoe',
				      -exportselection => 0,
				      -height => 6,
				      -selectmode => 'extended'
				      )->pack(-side => 'left',
					      -fill => 'both',
					      -expand => 1,
					      -padx => 0,
					      -pady => 0);
    my @x_vals = Eval::get_graph_x_types();
    foreach my $x (@x_vals){	
	$xlb->insert('end',$x);
    }
    $xlb->selectionSet(0);
    my @x_levels = Eval::get_graph_x_levels();
    my %y_types = Eval::get_graph_y_types();
    my @stat_lbs;
    my @type_map;
    foreach my $level (@x_levels){
	$x_level_lb->insert('end',$level);
    }	
    $x_level_lb->selectionSet(0);
    my $show_frame = $graph_graphs_frame->Frame()->pack(-side => 'bottom',
							-fill => 'x');
    $show_frame->Label(-text => "Selected Graphs:")->pack(-side => 'top');
    
    my $glb = $show_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				    -exportselection => 0,
				    -selectmode => 'extended'
				    )->pack(-side => 'top',
					    -fill => 'both',
					    -expand => 1,
					    -padx => 20,
					    -pady => 10);

    my $pslb = $select_pred_frame->Scrolled("Listbox", -scrollbars => 'osoe',
					    -exportselection => 0,
					    -height => 5,
					    -selectmode => 'extended'
					    )->pack(-side => 'top',
						    -fill => 'both',
						    -expand => 1,
						    -padx => 20,
						    -pady => 10);
    my $gslb = $select_graph_type_frame->Scrolled("Listbox",-scrollbars => 'osoe',
						  -exportselection => 0,
						  -height => 5,
						  -selectmode => 'single'
						  )->pack(-side => 'top',
							  -fill => 'both',
							  -expand => 1,
							  -anchor => 'w',
							  -padx => 20,
							  -pady => 10);
    my $level_frame = $select_graph_frame->Frame()->pack(-side => 'bottom',
							 -fill => 'both',
							 -expand => 0,
							 -anchor => 'w');
    my %sel_data;
    my %level_frames;


    #my %level_lbs;
    #my %type_lbs;
    #my %stat_lbs;
    foreach my $level (@x_levels){
	my %type_stat_frames;
	
	my $frame = $level_frame->Frame();
	$level_frames{$level} = $frame;
	
	my $level_lb = $frame->Scrolled("Listbox",-scrollbars => 'osoe',
					-exportselection => 0,
					-width => 10,
					-height => 5,
					-selectmode => 'single')->pack(-side => 'left',
								       -fill => 'y');
	my $below = 0;
	my $top_level;
	foreach my $l2 (@x_levels){
	    if($l2 eq $level){
		$below = 1;
	    }
	    if($below){
		unless(defined($top_level)){
		    $top_level = $l2;
		}
		$level_lb->insert('end',$l2);
		my $ts_frame = $frame->Frame();
		my $type_lb = $ts_frame->Scrolled("Listbox",-scrollbars => 'osoe',
						  -exportselection => 0,
						  -width  => 10,
						  -height => 5,
						  -selectmode => 'single'
						  )->pack(-side => 'left',
							  -fill => 'y',
							  -expand => 0);
		my @type_list = @{$y_types{$l2}{Type}};
		foreach my $type (@type_list){
		    $type_lb->insert('end',$type);
		}
		$type_lb->selectionSet(0);
		$sel_data{$level}{$l2}{type} = $type_list[0];
		$type_lb->bind("<Button>",sub{
		    my ($lb) = @_;
		    my @sel = $lb->curselection(); 
		    my $new_type = $lb->get($sel[0]);
		    $sel_data{$level}{$l2}{type} = $new_type;
		});
		#$type_lbs{$l2} = $type_lb;
		my $stat_lb = $ts_frame->Scrolled("Listbox",-scrollbars => 'osoe',
						  -exportselection => 0,
						  -width  => 5,
						  -height => 5,
						  -selectmode => 'single'
						  )->pack(-side => 'left',
							  -fill => 'both',
							  -expand => 1);
		my @stat_list = @{$y_types{$l2}{Stat}};
		foreach my $stat (@stat_list){
		    $stat_lb->insert('end',$stat);
		}
		$stat_lb->selectionSet(0);
	        $sel_data{$level}{$l2}{stat} = $stat_list[0];
		$stat_lb->bind("<Button>",sub{
		    my ($lb) = @_;
		    my @sel = $lb->curselection(); 
		    my $new_stat = $lb->get($sel[0]);
		    $sel_data{$level}{$l2}{stat} = $new_stat;
		});
		#$stat_lbs{$level} = $stat_lb;
		$type_stat_frames{$l2} = $ts_frame;
	    }
	}
	$level_lb->selectionSet(0);
	$type_stat_frames{$top_level}->pack(-fill => 'both',
					    -expand => 1);
	#$level_lbs{$level} = $level_lb;
	$sel_data{$level}{level} = $top_level;
	#my $cur_level = $x_levels[0];
	#$cur_levels{$level} = \$cur_level;
	$level_lb->bind("<Button>",sub{
	    my ($lb) = @_;
	    my @sel = $lb->curselection(); 
	    my $new_level = $lb->get($sel[0]);
	    my $cur_level = $sel_data{$level}{level};
	    #print "$cur_level = > $new_level\n";
	    if($cur_level ne $new_level){
		$type_stat_frames{$cur_level}->packForget();
		$type_stat_frames{$new_level}->pack(-fill => 'both',
						    -expand => 1);
		$sel_data{$level}{level} = $new_level;
	    }
	});
    }
    my $cur_level = "Gene";
    #EVAN need to not init to gene but instead init on gslb selectionSet
    my @sel = $gslb->curselection();
    if($#sel >= 0){
	$cur_level = $gslb->get($sel[0]);
	if($cur_level =~ /^(\S+)\:\:/){
	    $cur_level = $1;
	}
	else{
	    #ERROR
	}
    }
    $gslb->bind("<Button>",sub{
	my ($lb) = @_;
	my @sel = $lb->curselection(); 
	my $new_level = $lb->get($sel[0]);
	if($new_level =~ /^(\S+)\:\:/){
	    $new_level = $1;
	}
	if($cur_level ne $new_level){
	    $level_frames{$cur_level}->packForget();
	    $level_frames{$new_level}->pack(-fill => 'both',
					 -expand => 1);
	    $cur_level = $new_level;
	}
    });
    $gslb->bind("<Configure>",sub{
	my ($lb) = @_;
	my @sel = $lb->curselection(); 
	my $new_level = $lb->get($sel[0]);
	if($new_level =~ /^(\S+)\:\:/){
	    $new_level = $1;
	}
	if($cur_level ne $new_level){
	    $level_frames{$cur_level}->packForget();
	    $level_frames{$new_level}->pack(-fill => 'both',
					 -expand => 1);
	    $cur_level = $new_level;
	}
    });
    $level_frames{$cur_level}->pack(-fill => 'both',
				    -expand => 1);
    my $graphs;
    $ann_pred_button_frame->Button(-text => "Select Graphs ->",
				   -command => 
				   sub{
				       my @a_sel = $alb->curselection();
				       my @p_sel = $plb->curselection();
				       unless(($#a_sel >= 0) && ($#p_sel >= 0)){
					   error_func("Must select at least one ".
						      "annotation and prediction.");
					   return;
				       }
				       $ann_pred_frame->packForget();
				       $pslb->delete(0,'end');
				       foreach my $pred ($plb->curselection){
					   $pslb->insert('end',$plb->get($pred));
				       }
				       $pslb->selectionSet(0);
				       $graph_frame->pack(-side => 'top',
							  -fill => 'both',
							  -expand => 1); 
				   })->pack(-side => 'right');
    
    $graph_button_frame->Button(-text => "<- Choose Comparison",
				-command => 
				sub{
				    $graph_frame->packForget();
				    $ann_pred_frame->pack(-side => 'top',
							  -fill => 'both',
							  -expand => 1); 
				})->pack(-side => 'left');
    $graph_button_frame->Button(-text => "Add",
				-command => 
				sub{
				    my @split_sel = $xlb->curselection(); 
				    my @level_sel = $x_level_lb->curselection();
				    foreach my $x (@split_sel){
					my $split = $xlb->get($x);
					foreach my $l (@level_sel){
					    my $level = $x_level_lb->get($l);
					    my $graph = "$level"."::$split";
					    my $skip = 0;
					    foreach my $g ($glb->get(0,'end')){
						if($g eq $graph){
						    $skip = 1;
						}
					    }
					    unless($skip){
						$glb->insert('end',$graph);
					    }
					}
				    }
				})->pack(-side => 'left');
    
    $graph_button_frame->Button(-text => "Remove",
				-command => 
				sub{
				    my @sel = $glb->curselection();
				    while($#sel >= 0){
					$glb->delete($sel[0]);
					@sel = $glb->curselection();
				    }
				})->pack(-side => 'left');
    
    $graph_button_frame->Button(-text => "Create Graphs ->",
				-command => 
				sub{
				    my @sel = $glb->get(0,'end');
				    unless($#sel >= 0){
					error_func("Must add at least one graph type.");
					return;
				    }
				    $graphs = graph_run_func($alb,$plb,$glb,
							     \%Graph_Resolution);
				    $graph_frame->packForget();
				    $gslb->delete(0,'end');
				    my @graph_list = $glb->get(0,'end');
				    foreach my $graph (@graph_list){
					$gslb->insert('end',$graph);
				    }
				    $gslb->selectionSet(0);
				    $select_frame->pack(-side => 'top',
							-fill => 'both',
							-expand => 1); 
				}
				)->pack(-side => 'right');
    
    $select_button_frame->Button(-text => "<- Choose Graphs",
				 -command => 
				 sub{				    
				     $select_frame->packForget();
				     $graph_frame->pack(-side => 'top',
							-fill => 'both',
							-expand => 1); 
				 }
				 )->pack(-side => 'left');
    
    $select_button_frame->Button(-text => "Save",
				 -command => 
				 sub{
				     my $level = $sel_data{$cur_level}{level};
				     my $type = $sel_data{$cur_level}{$level}{type};
				     my $stat = $sel_data{$cur_level}{$level}{stat};
				     save_graph_func($alb,$pslb,$gslb,$graphs,
						     $level,$type,$stat)}
				 )->pack(-side => 'right');
    $select_button_frame->Button(-text => "View",
				 -command => 
				 sub{
				     my $level = $sel_data{$cur_level}{level};
				     my $type = $sel_data{$cur_level}{$level}{type};
				     my $stat = $sel_data{$cur_level}{$level}{stat};
				     display_graph_func($alb,$pslb,$gslb,$graphs,
						     $level,$type,$stat)}
				 )->pack(-side => 'right');
}

sub get_default_graph_resolution{
    return (uniform => {count => 10,
			min => 0,
			size => .1});
}

sub graph_run_func{
    my ($alb,$plb,$glb,$resolution) = @_;
    my @ann_index = $alb->curselection();
    my @pred_index = $plb->curselection();
    my $ann = $Gtf_Objs[$ann_index[0]];
    my $ann_name = $alb->get($ann_index[0]);
    my @preds;
    my @pred_names;
    foreach my $p_index (@pred_index){
	push @preds, $Gtf_Objs[$p_index];
	push @pred_names, $Obj_Names[$p_index];
    }
    my @graph_info;
    foreach my $graph ($glb->get(0,'end')){
	if($graph =~ /^(.+)::(.+)$/){	    
	    push @graph_info, {level => $1,
			       split => $2};
	}
	else{
	    error_func("Bad value in graph listbox in graph_run_func()",1);
	}
    }
    my @graphs = 
      Eval::make_graphs($ann,\@preds,\@graph_info,$resolution,$Verbose);
    return \@graphs;
}

sub save_graph_func{
    my ($alb,$pslb,$glb,$graphs,$level,$type,$stat) = @_;
    my @asel = $alb->curselection();
    my @psel = $pslb->curselection();
    my @gsel = $glb->curselection();
    my @pred_names;
    foreach my $index (@psel){
	push @pred_names, $pslb->get($index);
    }
    my $ann_name = $alb->get($asel[0]);
    my $temp = $glb->get($gsel[0]);
    unless($temp =~ /^(\S+)\:\:(\S+)$/){
	return;
    }
    my $split_type = $2;
    my $split_level = $1;
    my $filename = 
	$mw->getSaveFile(-defaultextension => ".graph",
			 -filetypes        =>
			 [['Graph Files', '*.graph'],
			  ['All Files','*']],
			 -initialdir       => $Cwd,
			 -initialfile      => "",
			 -title            => "Save Graph File");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
    }
    else{
	message_func("No file selected.  Graph not saved\n");
	return;
    }
    unless(open(OUT,">$filename")){
	error_func("Could not open \"$filename\" for write.  ".
		   "File not saved.");
	return;
    }
    print OUT "$level $type $stat vs $split_type\n";
    foreach my $pred (@pred_names){
	print OUT "\t$pred";
    }
    print OUT "\n";
    for(my $bin = 0;$bin <= $#{$$graphs[0]{$split_type}{$split_level}};$bin++){
	print OUT "".$$graphs[0]{$split_type}{$split_level}[$bin]{min}." - ".
	    $$graphs[0]{$split_type}{$split_level}[$bin]{max};
	foreach my $pred (@psel){
	    my $format = $$graphs[$pred]{$split_type}{$split_level}[$bin];
	    print OUT "\t".$$format{data}{$level}{$type}{$stat};
	}
	print OUT "\n";
    }
    close(OUT);
}

sub display_graph_func{
    my ($alb,$pslb,$glb,$graphs,$level,$type,$stat) = @_;
    my @asel = $alb->curselection();
    my @psel = $pslb->curselection();
    my @gsel = $glb->curselection();
    my @pred_names;
    foreach my $index (@psel){
	push @pred_names, $pslb->get($index);
    }
    my $ann_name = $alb->get($asel[0]);
    my $temp = $glb->get($gsel[0]);
    unless($temp =~ /^(\S+)\:\:(\S+)$/){
	return;
    }
    my $split_type = $2;
    my $split_level = $1;
    my @points;    
    for(my $bin = 0;$bin <= $#{$$graphs[0]{$split_type}{$split_level}};$bin++){
	my $min = $$graphs[0]{$split_type}{$split_level}[$bin]{min};
	my $max = $$graphs[0]{$split_type}{$split_level}[$bin]{max};
	foreach my $pred (@psel){
	    my $format = $$graphs[$pred]{$split_type}{$split_level}[$bin];
	    push @points, {min => $min,
			   max => $max,
			   count => $$format{data}{$level}{$type}{$stat}};
	}
    }
    gnuplot_bar_bin(\@points,"$level $stat $type vs $split_type","$split_type",
		    "$level $stat $type");
}

sub gnuplot_bar_bin{
    my ($graph,$title,$x_label,$y_label) = @_;
    my $tmp_file = get_temp_file();
    unless(open(TMP,">$tmp_file")){
	error_func("Couldn't write to temp file, $tmp_file,".
		   "for gnuplot.\n");
	return;
    }
    $graph = [sort {$$a{min} <=> $$b{min}} @$graph];
    my $x_min = $$graph[0]{min};
    my $x_max = $$graph[$#$graph]{max};    
    if($x_max =~ /\+/){
	$x_max = $$graph[$#$graph]{min};    
	my $size = $$graph[0]{max} - $$graph[0]{min};
	$x_max += $size;
	$$graph[$#$graph]{max} = $x_max;
    }
    my $y_max = -1;
    foreach my $point (@$graph){
	print TMP "$$point{min}\t0\n";
	print TMP "$$point{min}\t$$point{count}\n";
	print TMP "$$point{max}\t$$point{count}\n";
	print TMP "$$point{max}\t0\n";
	if($$point{count} > $y_max){
	    $y_max = $$point{count};
	}
    }
    close(TMP);
    if($y_max == 0){
	$y_max = 1;
    }
    my $size = $x_max - $x_min;
    my $buffer = $size *.05;
    $x_max += $buffer;
    $x_min -= $buffer;
    $y_max *= 1.05;
    if(open(GNU,"| $GNUPLOT -persist")){
	print GNU "set title \"$title\"\n";
	print GNU "set xlabel \"$x_label\"\n";
	print GNU "set ylabel \"$y_label\"\n";
	print GNU "plot [$x_min to $x_max] [0 to $y_max] \"$tmp_file\" with lines\n";
	print GNU "exit\n";
	close(GNU);
	system("rm $tmp_file");
    }
    else{
	error_func("gnuplot not found at $GNUPLOT");
    }
}

#############################################################
#
#  Overlap Frame Functions
# 
#############################################################
# initializes the overlap stats frame
sub initialize_overlap_frame{
    my $main_frame = $Main_Frames[$OVERLAP_FRAME_NUM]->Frame(-borderwidth => 2,
							     -relief => 'ridge'
							     )->pack(-side => 'top',
								     -expand => 1,
								     -fill => 'both');
    my ($frame,$plb) = make_pred_frame($main_frame);
    my $button_frame = $Main_Frames[$OVERLAP_FRAME_NUM]->Frame(-borderwidth => 2,
							       -relief => 'ridge'
							       )->pack(-side => 'bottom',
								       -fill => 'x');
    my $select_frame = $main_frame->Frame()->pack(-side => 'top',
						  -fill => 'x');
    $select_frame->Label(-text => "Choose overlap type:")->pack(-side => 'top');
    my $slb = $select_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				      -exportselection => 0,
				      -selectmode => 'single')->pack(-side => 'top',
								     -fill => 'x',
								     -expand => 0,
								     -padx => 20,
								     -pady => 10);
    my @overlap_list = Eval::get_overlap_mode_list();
    foreach my $overlap_type (@overlap_list){
	my $display_type = $overlap_type;
	$display_type =~ s/_/ /g;
	$slb->insert('end',$display_type);
    }
    $slb->selectionSet(0);
    $button_frame->Button(-text => "Get Overlap",
			  -command => sub{overlap_stats_run_func($plb,$slb)}
			  )->pack(-side => 'right');
}

sub overlap_stats_run_func{
    my ($plb,$slb) = @_;
    my @pred_index = $plb->curselection();
    my @preds;
    my @pred_names;
    foreach my $p_index (@pred_index){
	push @preds, $Gtf_Objs[$p_index];
	push @pred_names, $Obj_Names[$p_index];
    }
    my @mode_index = $slb->curselection();
    my $mode = $slb->get($mode_index[0]);
    $mode =~ s/ /_/g;
    unless($#preds >= 0){
	error_func("Must select at least one prediction set.");
	return;
    }
    my %data = Eval::get_overlap_statistics(\@preds,$mode,$Verbose);
    display_overlap_stats_func(\%data,\@pred_names,$mode);
}

sub display_overlap_stats_func{
    my ($data,$pred_names,$overlap_type) = @_;
    $overlap_type =~ s/_/ /g;
    my $window = $mw->Toplevel(-title => "Overlap Statistics");
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = $mw->height;
    my $width = $mw->width;
    my $report = "";
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $main = $window->Frame()->pack(-expand => 1,
				      -fill => 'both');
    
    my $main_pane = $main->Scrolled("Pane", 
				    -scrollbars => 'osoe',
				    -sticky => 'nswe')->pack(-fill => 'both',
							     -expand => 1);
    my $main_frame = $main_pane->Frame(-label => "Overlap Statistics",
				       -borderwidth => 2,
				       -relief => 'ridge')->pack(-side => 'top',
								 -expand => 1,
								 -fill => 'both');
    my $button_frame = $window->Frame(-borderwidth => 2,
				      -relief => 'ridge')->pack (-side =>'bottom',
								 -fill => 'x');
    $button_frame->Button(-text => 'Close',
			  -command => sub {$window->destroy}
			  )->pack(-side => 'right');
    $button_frame->Button(-text => 'Save',
			  -command => sub {
			      save_overlap_stats_output($data,$pred_names,
							$overlap_type)}
			  )->pack(-side => 'left');		
    my $overlap_frame = $main_frame->Frame(-borderwidth => 2,
					   -relief => 'ridge')->pack(-expand => 1,
								     -fill => 'both');
    my $type_frame = $overlap_frame->Frame(-borderwidth => 2,
					   -relief => 'ridge',
					   )->pack(-side => 'top',
						   -fill => 'x');
    $type_frame->Label(-text => "Overlap Type: $overlap_type",
			 -borderwidth => 2,
			 -relief => 'ridge')->pack(-side => 'top',
						   -fill => 'x');
    my $header_frame = $overlap_frame->Frame(-borderwidth => 2,
					     -relief => 'ridge',
					     )->pack(-side => 'top',
						     -fill => 'x');
    $header_frame->Label(-text => "Name (Count)",
			 -borderwidth => 2,
			 -relief => 'ridge')->pack(-side => 'top',
						   -fill => 'x');
    my $data_frame = $overlap_frame->Frame(-borderwidth => 2,
					   -relief => 'ridge'
					   )->pack(-side => 'top',
						   -fill => 'x');
    $data_frame->Label(-text => "Data",
		       -borderwidth => 2,
		       -relief => 'ridge')->pack(-side => 'top',
						 -fill => 'x');
    my %total = (all => 0);
    my @labels = Eval::get_overlap_labels($#$pred_names);
    foreach my $label (@labels){
	$total{$label} = 0;
    }    
    foreach my $group (keys %$data){
	$total{all} += $$data{$group}{total};
	foreach my $label (@labels){
	    $total{$label} += $$data{$group}{$label};
	}
    }
    for(my $i = 0;$i <= $#labels;$i++){
	$header_frame->Label(-text => 
			     "$labels[$i] -> $$pred_names[$i] ($total{$labels[$i]})"
			     )->pack(-side => 'top',
				     -anchor => 'w');
    }
    my $label_frame = $data_frame->Frame()->pack(-side => 'left');
    my $count_frame = $data_frame->Frame()->pack(-side => 'left');
    my @info_frames;
    for(my $i = 0;$i <= $#labels;$i++){
	my $frame = $data_frame->Frame()->pack(-side => 'left');
	push @info_frames, $frame;
    }
    my $total_frame = $data_frame->Frame()->pack(-side => 'left');    
    $label_frame->Label(-text => "Name")->pack(-side => 'top');
    $count_frame->Label(-text => "Count")->pack(-side => 'top');
    for(my $i = 0;$i <= $#labels;$i++){
	$info_frames[$i]->Label(-text => "$labels[$i]%"
				)->pack(-side => 'top');
    }
    $total_frame->Label(-text => "Total%")->pack(-side => 'top');
    foreach my $group (sort {length($a) <=> length($b) ||
				 $a cmp $b} (keys %$data)){
	$label_frame->Label(-text => $group)->pack(-side => 'top');
	$count_frame->Label(-text => $$data{$group}{total})->pack(-side => 'top');
	for(my $i = 0;$i <= $#labels;$i++){
	    my $val = sprintf $Precision, (100*$$data{$group}{$labels[$i]}/
					   $total{$labels[$i]});
	    $info_frames[$i]->Label(-text => "$val%")->pack(-side => 'top');
	}
	my $val = sprintf $Precision, (100*$$data{$group}{total}/$total{all});
	$total_frame->Label(-text => "$val%")->pack(-side => 'top');
    }
}

sub save_overlap_stats_output{
    my ($data,$pred_names,$overlap_type) = @_;
    my $report = get_overlap_stats_text($data,$pred_names,$overlap_type);
    my $filename = $mw->getSaveFile(-defaultextension => ".overlap",
				    -filetypes        =>
				    [['Overlap Output','.overlap'],
				     ['All Files','*']],
				    -initialdir       => $Cwd,
				    -initialfile      => "",
				    -title            => "Save Overlap Output As");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
	if(open(OUT,">$filename")){
	    print OUT $report;
	    message_func("File \"$filename\" saved.");
	    close(OUT);
	}
	else{
	    error_func("Could not open \"$filename\" for write.");
	}
    }
    else{
	error_func("No file selected.  File not saved");
    }
    
}

sub get_overlap_stats_text{
    my ($data,$pred_names,$overlap_type) = @_;
    my $report = "Overlap Type: $overlap_type\n";
    my %total = (all => 0);
    my @labels = Eval::get_overlap_labels($#$pred_names);
    foreach my $label (@labels){
	$total{$label} = 0;
    }    
    foreach my $group (keys %$data){
	$total{all} += $$data{$group}{total};
	foreach my $label (@labels){
	    $total{$label} += $$data{$group}{$label};
	}
    }
    for(my $i = 0;$i <= $#labels;$i++){
	$report .= "$labels[$i] -> $$pred_names[$i] ($total{$labels[$i]})\n";
    }
    $report .= "\nname\tcount";
    for(my $i = 0;$i <= $#labels;$i++){
	$report .= "\t$labels[$i]%";
    }
    $report .= "\ttotal\n";
    foreach my $group (sort {length($a) <=> length($b) ||
				 $a cmp $b} (keys %$data)){
	$report .= $group ."\t".$$data{$group}{total};
	for(my $i = 0;$i <= $#labels;$i++){
	    $report .= sprintf "\t".$Precision, (100*$$data{$group}{$labels[$i]}/
						 $total{$labels[$i]});
	    $report .= "%";
	}
	$report .= sprintf "\t".$Precision, (100*$$data{$group}{total}/$total{all});
	$report .= "%\n";
    }
    return $report;
}

#############################################################
#
#  Distribution Frame Functions
# 
#############################################################
# initializes the distribution stats frame
sub initialize_distribution_frame{
    my $main_frame = $Main_Frames[$DIST_FRAME_NUM]->Frame()->pack(-expand => 1,
								  -fill => 'both');
    my $select_frame = $main_frame->Frame()->pack(-side => 'top',
						  -fill => 'both',
						  -expand => 1);
    my $display_frame = $main_frame->Frame();
    my ($pred_frame,$plb) = make_pred_frame($select_frame);
    $pred_frame->configure(-borderwidth => 2,
			   -relief => 'ridge');
    $pred_frame->pack(-side => 'top',
		      -fill => 'both',
		      -expand => 1);
    $pred_frame->Label(-text => "Select Distribution Type:")->pack();
    my $dlb = $pred_frame->Listbox(-exportselection => 0,
				   -selectmode => 'extended')->pack(-fill => 'both',
								    -padx => 20,
								    -pady => 10);
    my @distribution_names = Eval::get_distribution_type_list();
    my %distributions = Eval::get_distribution_type_hash();
    foreach my $dist (@distribution_names){
	$dist =~ s/_/ /g;
	$dlb->insert('end',$dist);
    }
    $dlb->selectionSet(0);
    $display_frame->Label(-text => 'Choose Distribution')->pack(-side => 'top',
								-fill => 'y');
    my $display_lb = $display_frame->Listbox(-exportselection => 0,
					     -selectmode => 'single'
					     )->pack(-side => 'top',
						     -fill => 'both',
						     -padx => 20,
						     -pady => 10);   
    my $cum_frame = $display_frame->Frame()->pack(-side => 'top',
						  -fill => 'x',
						  -padx => 20);
    $cum_frame->Label(-text => "Cumulative Distribution:"
		      )->pack(-side => 'left');
    my $cum = 0;
    $cum_frame->Checkbutton(-variable => \$cum)->pack(-side => 'right');
    my $max_frame = $display_frame->Frame()->pack(-side => 'top',
						  -fill => 'x',
						  -padx => 20);
    $max_frame->Label(-text => "Select upper bound of distribution:"
		      )->pack(-side => 'left');
    my $max = 1000;
    $max_frame->Entry(-textvariable => \$max)->pack(-side => 'right');
    my $res_frame = $display_frame->Frame()->pack(-side => 'top',
						  -fill => 'x',
						  -padx => 20);
    $res_frame->Label(-text => "Select bin size:"
		      )->pack(-side => 'left');
    my $res = 50;
    $res_frame->Entry(-textvariable => \$res)->pack(-side => 'right');
    my $select_button_frame = $select_frame->Frame(-borderwidth => 2,
						   -relief => 'ridge'
						   )->pack(-side => 'bottom',
							   -fill => 'x');
    my @data;
    $select_button_frame->Button(-text => "Get Distribution->",
				 -command => sub{
				     my @pred_index = $plb->curselection();
				     my @preds;
				     my @pred_names;
				     foreach my $p_index (@pred_index){
					 push @preds, $Gtf_Objs[$p_index];
					 push @pred_names, $Obj_Names[$p_index];
				     }
				     foreach my $type (keys %distributions){
					 $distributions{$type} = 0;
				     }
				     my @dist_index = $dlb->curselection();
				     foreach my $d_index (@dist_index){
					 my $dist_type = $dlb->get($d_index);
					 $dist_type =~ s/ /_/g;
					 $distributions{$dist_type} = 1;
				     }
				     unless($#preds >= 0){
					 error_func("Must select at least one ".
						    "prediction set.");
					 return;
				     }
				     @data = 
				       Eval::get_distribution(\@preds,\%distributions);
				     $display_lb->delete(0,'end');
				     for(my $i = 0;$i <= $#data;$i++){
					 foreach my $dist (keys %{$data[$i]}){
					     $dist =~ s/_/ /g;
					     $display_lb->insert('end',
								 "$pred_names[$i]: $dist"
								 );
					 }
				     }
				     $display_lb->selectionSet(0);
				     $select_frame->packForget();
				     $display_frame->pack(-side => 'top',
							  -fill => 'both',
							  -expand => 1); 
				 }
				 )->pack(-side => 'right');
    my $display_button_frame = $display_frame->Frame(-borderwidth => 2,
						     -relief => 'ridge'
						     )->pack(-side => 'bottom',
							     -fill => 'x');
    $display_button_frame->Button(-text => "View",
				  -command => sub{
				      my @dist_index = $display_lb->curselection();
				      my $dist = $display_lb->get($dist_index[0]);
				      $dist =~ /^(.+): (.+)$/;
				      my $gtf_name = $1;
				      $dist = $2;				      
				      my $dist_key = $dist;
				      $dist_key =~ s/ /_/g;
				      my $i;
				      for($i = 0;$i <= $#Obj_Names;$i++){
					  if($Obj_Names[$i] eq $gtf_name){
					      last;
					  }
				      }
				      if($i > $#Obj_Names){
					  error_func("Internal Error: Can't locate GTF");
				      }
				      unless(($max =~ /\d+.\d+/) ||
					     ($max =~ /\d+/)){
					  error_func("Bad \"Max\" field.  Must be a ".
						     "number.");
					  return;
				      }
				      unless(($res =~ /\d+.\d+/) ||
					     ($res =~ /\d+/)){
					  error_func("Bad \"bin size\" field.  Must be ".
						     "a number.");
					  return;
				      }
				      if($cum){
					  $dist .= " Cumulative";
				      }
				      graph_distribution_data($gtf_name,$dist,
							     bin_distribution_data
							     ($data[$i]{$dist_key},$res,
							      $max,$cum));
				  }
				  )->pack(-side => 'right');
    $display_button_frame->Button(-text => "Save",
				  -command => sub{
				      my @dist_index = $display_lb->curselection();
				      my $dist = $display_lb->get($dist_index[0]);
				      $dist =~ /^(.+): (.+)$/;
				      my $gtf_name = $1;
				      $dist = $2;				      
				      my $dist_key = $dist;
				      $dist_key =~ s/ /_/g;
				      my $i;
				      for($i = 0;$i <= $#Obj_Names;$i++){
					  if($Obj_Names[$i] eq $gtf_name){
					      last;
					  }
				      }
				      if($i > $#Obj_Names){
					  error_func("Internal Error: Can't locate GTF");
				      }
				      unless(($max =~ /^\d+\.\d+$/) ||
					     ($max =~ /^\d+$/)){
					  error_func("Bad \"Max\" field.  Must be a ".
						     "number.");
					  return;
				      }
				      unless(($res =~ /^\d+\.\d+$/) ||
					     ($res =~ /^\d+$/)){
					  error_func("Bad \"bin size\" field.  Must be ".
						     "a number.");
					  return;
				      }
				      if($cum){
					  $dist .= " Cumulative";
				      }
				      save_distribution_data($gtf_name,$dist,
							     bin_distribution_data
							     ($data[$i]{$dist_key},$res,
							      $max,$cum));
				  }
				  )->pack(-side => 'right');
    $display_button_frame->Button(-text => "<-Back",
				  -command => sub{
				      $display_frame->packForget();
				      $select_frame->pack(-side => 'top',
							   -fill => 'both',
							   -expand => 1); 	
				  }
				  )->pack(-side => 'left');
    
    
}

sub save_distribution_data{
    my ($gtf_name,$dist_name,$data) = @_;
    my $dist_key = $dist_name;
    $dist_key =~ s/ /_/g;
    my $filename = 
	$mw->getSaveFile(-defaultextension => ".dist",
			 -filetypes        => [['Distribution Files',['.dist']],
					       ['All Files','*']],
			 -initialdir       => $Cwd,
			 -initialfile      => "",
			 -title            => "Save File");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
	unless(open(OUT,">$filename")){
	    error_func("Couldn't open $filename for write.  File not saved.");	
	    return;
	}
	print OUT "$dist_name\n";
	my $real_max = 0;
	foreach my $point (@$data){
	    if($$point{max} eq "+"){
		print OUT ">".$$point{min}."\t".$$point{count}."\n";
	    }
	    else{
		print OUT $$point{min}."-".$$point{max}."\t".$$point{count}."\n";
	    }
	}
	if($Verbose){
	    print STDERR "$filename saved.\n";
	}
	close(OUT);
    }
}

sub graph_distribution_data{
    my ($gtf_name,$dist_name,$data) = @_;
    my $dist_key = $dist_name;
    $dist_key =~ s/ /_/g;
    my $title = "$gtf_name - $dist_name Distribution";
    my $y_label = "Count";
    if($dist_name =~ /Gene/){
	$y_label = "Gene Count";
    }
    elsif($dist_name =~ /Transcript/){
	$y_label = "Transcipt Count";
    }
    elsif($dist_name =~ /Exon/){
	$y_label = "Exon Count";
    }    
    my $x_label = "$dist_name";
    gnuplot_bar_bin($data,$title,$x_label,$y_label);
    return;
}

sub bin_distribution_data{
    my ($data,$res,$max,$cum) = @_;
    my @vals = keys %$data;
    @vals = sort {$a <=> $b} @vals;
    my $bin_max = $res;
    my $bin_min = 0;
    my @output;
    my $count = 0;
    foreach my $val (@vals){
	while($val > $bin_max){
	    push @output, {min => $bin_min,
			   max => $bin_max,
			   count => $count};
	    unless($cum){
		$count = 0;
	    }
	    $bin_min = $bin_max;
	    $bin_max += $res;
	    if(($bin_max > $max) && ($bin_min < $max)){
		$bin_max = $max;
	    }
	    elsif($bin_min == $max){
		$bin_max = $vals[$#vals];
	    }
	}
	$count += $$data{$val};
    }
    while($bin_min < $max){
	push @output, {min => $bin_min,
		       max => $bin_max,
		       count => $count};
	unless($cum){
	    $count = 0;
	}
	$bin_min = $bin_max;
	$bin_max += $res;
	if(($bin_max > $max) && ($bin_min < $max)){
	    $bin_max = $max;
	}
    }
    push @output, {min => $bin_min,
		   max => "+",
		   count => $count};
    return \@output;
}

#############################################################
#
#  General Functions
# 
#############################################################
sub make_ann_pred_frame{
    my ($frame) = @_;
    my $ann_pred_frame = $frame->Frame();
    my $ann_frame = $ann_pred_frame->Frame()->pack(-side => 'top',
						   -fill => 'x',
						   -expand => 0);
    $ann_frame->Label(-text => "Select Annotation:    ")->pack(-side => 'top');
    my $alb = $ann_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				   -exportselection => 0,
				   -selectmode => 'single')->pack(-side => 'top',
								  -fill => 'x',
								  -expand => 0,
								  -padx => 20,
								  -pady => 10);
    push @Ann_Gtf_Lbs, $alb;
    $ann_frame->Label(-text => " ")->pack(-side => "top");
    my ($pred_frame,$plb) =  make_pred_frame($ann_pred_frame);
    $pred_frame->pack(-side => 'top',
		      -fill => 'both',
		      -expand => 1);
    return ($ann_pred_frame,$alb,$plb);
}

sub make_pred_frame{
    my ($frame) = @_;
    my $pred_frame = $frame->Frame()->pack(-side => 'top',
					   -fill => 'both',
					   -expand => 1);
    $pred_frame->Label(-text => "Select Predictions:    ")->pack(-side => 'top');
    my $plb = $pred_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				    -exportselection => 0,
				    -selectmode => 'extended')->pack(-side => 'top',
								     -fill => 'both',
								     -expand => 1,
								     -padx => 20,
								     -pady => 10);
    $pred_frame->Label(-text => " ")->pack(-side => "top");
    push @Pred_Gtf_Lbs, $plb;
    return ($pred_frame,$plb);
}

sub adjust_data_precision{
    my ($data) = @_;
    my %order = Eval::get_list_struct();
    for(my $i = 0; $i <= $#$data;$i++){
	foreach my $level (@{$order{Levels}}){
	    foreach my $type (@{$order{$level}{Type}}){
		foreach my $stat (@{$order{$level}{Stat}}){
		    $$data[$i]{$level}{$type}{$stat} = 
			sprintf $Precision, $$data[$i]{$level}{$type}{$stat};
		    if(($stat =~ /Sensitivity/) ||
		       ($stat =~ /Specificity/)){
			$$data[$i]{$level}{$type}{$stat} .= "%";
		    }
		    else{
			$$data[$i]{$level}{$type}{$stat} .= " ";
		    }
		}
	    }
	}
    }    
}

# Display message to User
sub message_func{
    my ($msg) = @_;
    if($Verbose){
	print "$msg\n";
    }
    $mw->messageBox(-message => $msg, 
		    -type => 'OK', 
		    -default => 'OK'); 
}

# Display Error to User
sub error_func{
    my ($error_msg,$fatal) = @_;
    unless(defined($fatal)){
	$fatal = 0;
    }
    if($fatal){
	$error_msg .= "\nExiting.";
    }
    $mw->messageBox(-icon => 'questhead', 
		    -message => $error_msg, 
		    -type => 'OK', 
		    -default => 'OK'); 
    if($fatal){
	exit($fatal);
    }
}

sub get_temp_file{
    $Next_Temp++;
    my $pid = $$;
    return "/tmp/$USER.$pid.$Next_Temp.temp";
}

#############################################################
#
#  Menu Functions
# 
#############################################################
sub open_func{
    my $filename = 
	$mw->getOpenFile(-defaultextension => ".gtf",
			 -filetypes        =>
			 [['GTF Files',['.gtf','.gff']],
			  ['List Files','.list'],
			  ['All Files','*']],
			 -initialdir       => $Cwd,
			 -initialfile      => "",
			 -title            => "Open File");
    if(defined($filename)){
	if($filename =~ /^(.+)\/\S+$/){
	    $Cwd = $1;
	}
	load_func($filename);
    }
}

sub load_func{
    my ($filename) = @_;
    if(-e $filename){
	my $gtf;
	if(($filename =~ /^.*\.list/) || 
	   ($List_Mode && ($filename !~/^.*\.g[tf]f$/))){
	    open(LIST,$filename) 
		or error_func("Unable to open $filename for input.",1);
	    while(my $file = <LIST>){
		if($file =~ /^\s+(\S+)\s+$/){
		    $file = $1;
		}
		elsif($file =~ /^\s+(\S+)$/){
		    $file = $1;
		}
		elsif($file =~ /^(\S+)\s+$/){
		    $file = $1;
		}
		if(($file =~ /\S/) &&
		   ($file !~ /^\#/)){
		    chomp $file;
		    my @info = split /\t/, $file;			
		    my $next = create_gtf_object($info[0],$info[1],$info[2]);
		    push @$gtf, $next;
		}
	    }
	    close(LIST);
	}
	else{
	    $gtf = [create_gtf_object($filename)];		
	}
	add_to_display($filename,$gtf);
    }
    else{
	error_func("File \"$filename\" does not exist.");
    }
}

sub create_gtf_object{
    #sequence autoloading has been commented out below
    my ($file,$seq,$conseq) = @_;
    my $info = {gtf_filename => $file};
    if($Really_Verbose){
	$info->{warning_fh} = \*STDERR;
    }	
    unless($No_Seq){
	my $base = $file;
	if($file =~ /^(.+)\.g.f/i){
	    $base = $1;
	}
	if(defined($conseq) && length($conseq) > 0){
	    $info->{conseq_filename} = $conseq;
	}  
	elsif(-e "$base.conseq"){
	    #$info->{conseq_filename} = "$base.conseq";
	}
	if(defined($seq) && length($seq) > 0){
	    $info->{seq_filename} = $seq;
	}  
	elsif(-e "$base.seq"){
	    #$info->{seq_filename} = "$base.seq";
	}
	elsif(-e "$base.fasta"){
	    #$info->{seq_filename} = "$base.fasta";
	}
	elsif(-e "$base.fa"){
	    #$info->{seq_filename} = "$base.fa";
	}
	elsif((-e "$base") && (($base =~ /\.seq$/) or 
			       ($base =~ /\.fasta$/) or
			       ($base =~ /\.fa$/))){
	    #$info->{seq_filename} = "$base";
	}
    }
    if($Verbose){
	print STDERR "Loading $file";
	if((defined($info->{seq_filename})) and 
	   (defined($info->{conseq_filename}))){
	    print STDERR " with conservation sequence and fasta file.\n";
	}
	elsif(defined($info->{seq_filename})){
	    print STDERR " with fasta file.\n";
	}
	elsif(defined($info->{conseq_filename})){
	    print STDERR " with conservation sequence.\n";
	}
	else{
	    print STDERR "\n";
	}
    }
    if($Quick_Load){
	$info->{no_check} = 1;
    }
    return GTF::new($info);
}

# Add newly loaded gtf to all gtf selection list boxes
sub add_to_display{
    my ($name,$gtf) = @_;
    foreach my $lb (@Pred_Gtf_Lbs){
	$lb->insert('end',$name);
	my @sel = $lb->curselection;
	if($#sel == -1){
	    $lb->selectionSet(0);
	}
    }
    foreach my $lb (@Ann_Gtf_Lbs){
	$lb->insert('end',$name);
	my @sel = $lb->curselection;
	if($#sel == -1){
	    $lb->selectionSet(0);
	}
    }
    push @Gtf_Objs, $gtf;
    push @Obj_Names, $name;
}


sub save_func{
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = int($mw->height/2);
    my $width = 247;
    my $window = $mw->Toplevel(-title => "Save Predictions");
    $window->resizable(0,1);
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $window_pane = $window->Scrolled("Pane", 
					-scrollbars => 'oe',
					-sticky => 'nswe',
					);
    
    my $display_frame = $window_pane->Frame()->pack(-fill => 'both',
						    -expand => 1);
    
    $window_pane->pack(-expand => 1, -fill => 'both');
    my $button_frame = $window->Frame(-relief => "ridge",
				      -borderwidth => 2)->pack(-side => 'bottom', 
							       -fill => 'x');
    my $pred_frame = $display_frame->Frame()->pack(-side => 'top',
						   -fill => 'both',
						   -expand => 1);
    $pred_frame->Label(-text => "Select Predictions:    ")->pack(-side => 'top');
    my $plb = $pred_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				    -exportselection => 0,
				    -selectmode => 'extended')->pack(-side => 'top',
								     -fill => 'both',
								     -expand => 1,
								     -padx => 20,
								     -pady => 10);
    $pred_frame->Label(-text => " ")->pack(-side => "top");
    for(my $i = 0;$i <= $#Obj_Names;$i++){
	$plb->insert('end',$Obj_Names[$i]);
    }
    $button_frame->Button(-text => "Save",
			  -command => sub {
			      my @p_index = $plb->curselection();
			      foreach my $index (@p_index){
				  _save_pred_func($Gtf_Objs[$index]);
			      }
			  }
			  )->pack(-side => 'left');
    $button_frame->Button(-text => "Close",
			  -command => sub {$window->destroy})->pack(-side => 'right');
}

sub _save_pred_func{
    my ($object) = @_;
    if($#$object >= 1){
	my $filename = $mw->getSaveFile(-defaultextension => ".list",
					-filetypes        =>
					[['List','.list'],
					 ['All Files','*']],
					-initialdir       => $Cwd,
					-initialfile      => "",
					-title            => "Save List File As");
	if(defined($filename)){
	    if($filename =~ /^(.+)\/\S+$/){
		$Cwd = $1;
	    }
	    if(open(OUT,">$filename")){
		my $num = 0;
		foreach my $gtf (@$object){
		    $num++;
		    my $gtf_name = "$filename.$num.gtf";
		    unless(open(GTF,">$gtf_name")){
			error_func("File \"$gtf_name\" could not be opened for write.".
				   "  Files not saved\n");
			close(OUT);
			return;
		    }
		    $gtf->output_gtf_file(\*GTF);
		    close(GTF);
		    print OUT "$gtf_name\n";
		}
		message_func("File \"$filename\" saved.");
		close(OUT);
	    }
	    else{
		error_func("Could not open \"$filename\" for write.");
	    }
	}
	else{
	    error_func("No file selected.  File not saved");
	}
    }
    else{
	my $filename = $mw->getSaveFile(-defaultextension => ".gtf",
					-filetypes        =>
					[['GTF','.gtf'],
					 ['All Files','*']],
					-initialdir       => $Cwd,
					-initialfile      => "",
					-title            => "Save GTF File As");
	if(defined($filename)){
	    if($filename =~ /^(.+)\/\S+$/){
		$Cwd = $1;
	    }
	    if(open(OUT,">$filename")){
		$$object[0]->output_gtf_file(\*OUT);
		message_func("File \"$filename\" saved.");
		close(OUT);
	    }
	    else{
		error_func("Could not open \"$filename\" for write.");
	    }
	}
	else{
	    error_func("No file selected.  File not saved");
	}
    }
}

sub remove_func{
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = int($mw->height/2);
    my $width = 247;
    my $window = $mw->Toplevel(-title => "Remove Predictions");
    $window->resizable(0,1);
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $window_pane = $window->Scrolled("Pane", 
					-scrollbars => 'oe',
					-sticky => 'nswe',
					);
    
    my $display_frame = $window_pane->Frame()->pack(-fill => 'both',
						    -expand => 1);
    
    $window_pane->pack(-expand => 1, -fill => 'both');
    my $button_frame = $window->Frame(-relief => "ridge",
				      -borderwidth => 2)->pack(-side => 'bottom', 
							       -fill => 'x');
    my $pred_frame = $display_frame->Frame()->pack(-side => 'top',
						   -fill => 'both',
						   -expand => 1);
    $pred_frame->Label(-text => "Select Predictions:    ")->pack(-side => 'top');
    my $plb = $pred_frame->Scrolled("Listbox", -scrollbars => 'osoe',
				    -exportselection => 0,
				    -selectmode => 'extended')->pack(-side => 'top',
								     -fill => 'both',
								     -expand => 1,
								     -padx => 20,
								     -pady => 10);
    $pred_frame->Label(-text => " ")->pack(-side => "top");
    for(my $i = 0;$i <= $#Obj_Names;$i++){
	$plb->insert('end',$Obj_Names[$i]);
    }
    $button_frame->Button(-text => "Remove",
			  -command => sub {
			      _remove_preds_func($plb);
			      my @p_index = $plb->curselection;
			      while($#p_index >= 0){
				  $plb->delete($p_index[0]);
				  @p_index = $plb->curselection();
			      }
			  }
			  )->pack(-side => 'left');
    $button_frame->Button(-text => "Close",
			  -command => sub {$window->destroy})->pack(-side => 'right');
}

sub _remove_preds_func{
    my ($plb) = @_;
    my @p_index = $plb->curselection();
    my @rem_preds;
    foreach my $index (@p_index){
	$rem_preds[$index] = 1;
    }
    my @map;
    my $next = 0;
    for(my $i = 0;$i <= $#Gtf_Objs;$i++){
	if($rem_preds[$i]){
	    $map[$i] = 'X';
	}
	else{
	    $map[$i] = $next;
	    $next++;
	}
    }
    my @p_old_sel;
    foreach my $lb (@Pred_Gtf_Lbs){
	my @sel = $lb->curselection;
	push @p_old_sel, \@sel;
	$lb->delete(0,'end');
    }
    my @a_old_sel;
    foreach my $lb (@Ann_Gtf_Lbs){
	my @sel = $lb->curselection;
	push @a_old_sel, $sel[0];
	$lb->delete(0,'end');
    }
    my @new_gtfs;
    my @new_names;
    for(my $i = 0;$i <= $#Gtf_Objs;$i++){
	unless(defined($rem_preds[$i])){
	    push @new_gtfs, $Gtf_Objs[$i];
	    push @new_names, $Obj_Names[$i];
	    foreach my $lb (@Pred_Gtf_Lbs){
		$lb->insert('end',$Obj_Names[$i]);
		#$lb->selectionSet('end');
	    }
	    foreach my $lb (@Ann_Gtf_Lbs){
		$lb->insert('end',$Obj_Names[$i]);
	    }
	}
    }
    for(my $i = 0;$i <= $#Pred_Gtf_Lbs;$i++){
	foreach my $pred (@{$p_old_sel[$i]}){
	    unless($map[$pred] eq 'X'){
		$Pred_Gtf_Lbs[$i]->selectionSet($map[$pred]);
	    }
	}
	my @sel = $Pred_Gtf_Lbs[$i]->curselection;
	if($#sel == -1){
	    $Pred_Gtf_Lbs[$i]->selectionSet(0);
	}
    }
    for(my $i = 0;$i <= $#Ann_Gtf_Lbs;$i++){
	my $old_sel = $a_old_sel[$i];
	if($map[$old_sel] ne 'X'){
	    $Ann_Gtf_Lbs[$i]->selectionSet($map[$old_sel]);
	}
	else{
	    my @sel = $Ann_Gtf_Lbs[$i]->curselection;
	    if($#sel == -1){
		$Ann_Gtf_Lbs[$i]->selectionSet(0);
	    }
	}
    }
    @Gtf_Objs = @new_gtfs;
    @Obj_Names = @new_names;
}

sub help_func{
    if($Current_Frame == $EVAL_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening Eval screen help page.\n";
	}
	system("$GV $HELP_PATH/eval-help.ps &");
    }
    elsif($Current_Frame == $STATS_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening GenStats screen help page.\n";
	}
	system("$GV $HELP_PATH/stats-help.ps &");	
    }
    elsif($Current_Frame == $FILTER_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening Filter screen help page.\n";
	}
	system("$GV $HELP_PATH/filter-help.ps &");
    }
    elsif($Current_Frame == $GRAPH_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening Graph screen help page.\n";
	}
	system("$GV $HELP_PATH/graph-help.ps &");
    }
    elsif($Current_Frame == $OVERLAP_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening Overlap screen help page.\n";
	}
	system("$GV $HELP_PATH/overlap-help.ps &");
    }
    elsif($Current_Frame == $DIST_FRAME_NUM){
	if($Verbose){
	    print STDERR "Opening Dist screen help page.\n";
	}
	system("$GV $HELP_PATH/dist-help.ps &");
    }
}

sub about_func{
    my $window = $mw->Toplevel(-title => "About Eval");     
    my $title = $window->Frame(-label => "Eval",
			       -relief => 'ridge',
			       -borderwidth => 2)->pack(-side => "top",-fill => "x");
    my $write = $window->Frame(-relief => 'ridge',
			       -borderwidth => 2)->pack(-side => "top",-fill => "x");
    my $info = $window->Frame( -relief => 'ridge',
			       -borderwidth => 2)->pack(-side => "top",-fill => "x");
    $write->Label(-text => "Eval was written by Evan Keibler with the help of Michael Brent")->pack(-side => "top");
    $info->Label(-text => "For more information visit\nhttp://genome.cse.wustl.edu\n".
		   "then go to Resources then Software")->pack(-side => "top");

    
}

# Do any cleanup needed before termination here 
sub exit_func{
    exit;
}

#############################################################
#
#  Options Functions
# 
#############################################################
# load the options file
sub load_options_file{
    unless(-e $Options_File){
	if($Verbose){
	    print STDERR "Options file does not exist.  Creating.\n";
	}
	save_options_file();
    }
    unless(open(OPT,$Options_File)){
	error_func("Unable to open options preferences file for input.  ".
		   "Loading default options");
	return;
    }
    #load the file to %Display
    while(my $line = <OPT>){
	chomp $line;
	if($line eq "\@Display"){
	    while(my $display = <OPT>){
		chomp $display;
		if($display eq "\@Done"){
		    last;
		}
		$display =~ /^\t(.+)$/;
		$display = $1;
		my @info = split /\t/,$display;
		$Display{$info[0]}{$info[1]}{$info[2]} = $info[3];
	    }	    
	}
	elsif($line eq "\@Graph"){
	    while(my $display = <OPT>){
		chomp $display;
		if($display eq "\@Done"){
		    last;
		}
		$display =~ /^\t(.+)$/;
		$display = $1;
		my ($graph_type,$split_type,@info) = split /\t/,$display;
		if($split_type eq 'user'){
		    my @splits;
		    my $last_stop = $info[0];
		    for(my $i = 1;$i <= $#info;$i++){
			push @splits, {min => $last_stop,
				       max => $info[$i]};
			$last_stop = $info[$i];
		    }
		    $Graph_Resolution{$graph_type}{user} = \@splits;
		    if(defined($Graph_Resolution{$graph_type}{uniform})){
			delete($Graph_Resolution{$graph_type}{uniform});
		    }
		}
		elsif($split_type eq 'uniform'){
		    for(my $i = 0;$i <= 2;$i++){
			unless(defined($info[$i])){
			    $info[$i] = '';
			}
		    }
		    $Graph_Resolution{$graph_type}{uniform} = {min => $info[0],
							       size => $info[1],
							       count => $info[2]};
		    if(defined($Graph_Resolution{$graph_type}{user})){
			delete($Graph_Resolution{$graph_type}{user});
		    }
		}
		else{
		    error_func("Bad graph type, $graph_type, in options file.");
		}
	    }	    	    
	}
    }
    close(OPT);
}

# save the options file
sub save_options_file{
    unless(open(OPT,">$Options_File")){
	error_func("Unable to open options preferences file for output.  ".
		   "Options cannot be saved.");
	return;
    }
    #write %Display to disk
    print OPT "\@Display\n";
    foreach my $level (keys %Display){
	foreach my $catagory (keys %{$Display{$level}}){	    
	    foreach my $option (keys %{$Display{$level}{$catagory}}){
		print OPT "\t$level\t$catagory\t$option\t".
		    $Display{$level}{$catagory}{$option}."\n";
	    }
	}
    }
    print OPT "\@Done\n";
    print OPT "\@Graph\n";
    foreach my $type (keys %Graph_Resolution){
	print OPT "\t$type\t";
	foreach my $graph_type (keys %{$Graph_Resolution{$type}}){
	    if($graph_type eq 'uniform'){
		print OPT "uniform".
		    "\t$Graph_Resolution{$type}{$graph_type}{min}".
			"\t$Graph_Resolution{$type}{$graph_type}{size}\t".
			    "$Graph_Resolution{$type}{$graph_type}{count}\n";
	    }
	    elsif($graph_type eq 'user'){
		print OPT "user";
		print OPT "\t$Graph_Resolution{$type}{user}[0]{min}";
		foreach my $bin (@{$Graph_Resolution{$type}{user}}){
		    print OPT "\t$$bin{max}";
		}
		print OPT "\n";
	    }
	    else{
		error_func("Bad graph split type, $graph_type,".
			   "while saving options file.  This should never happen.");
	    }
	}
    }
    print OPT "\@Done\n";
    close(OPT);
}

#Load options screen
sub edit_options_func{
    my $t = "    ";
    my $start_x = $mw->x;
    my $start_y = $mw->y;
    my $height = $mw->height - 20;
    my $width = 247;
    my $window = $mw->Toplevel(-title => "Options");
    $window->resizable(0,1);
    $window->geometry("$width"."x$height+$start_x+$start_y");
    my $button_frame = $window->Frame(-relief => "ridge",
				      -borderwidth => 2)->pack(-side => 'bottom', 
							       -fill => 'x');
    my $select_frame = $window->Frame(-relief => 'ridge',
				      -borderwidth => 2)->pack(-side => 'top',
							       -fill => 'x'); 
    my @options_frames;
    my @options_buttons;
    my $eval_index = 0;
    my $graph_index = 1;
    my $current_frame = $eval_index;
    my $switch_func = sub {
	my ($new_frame) = @_;
	$options_frames[$current_frame]->packForget();
	$options_buttons[$current_frame]->configure(-relief => 'sunken',
						    -foreground => $INACTIVE_COLOR);
	$options_frames[$new_frame]->pack(-side => 'top',
					  -expand => 1,
					  -fill => 'both');
	$options_buttons[$new_frame]->configure(-relief => 'raised',
						-foreground => $ACTIVE_COLOR);
	$current_frame = $new_frame;
    };
    $options_buttons[$graph_index] = 
	$select_frame->Button(-text => '      Graph       ',
			      -relief => 'sunken',
			      -foreground => $INACTIVE_COLOR,
			      -command => sub{&$switch_func($graph_index)});
    $options_buttons[$eval_index] = 
	$select_frame->Button(-text => '    Eval Output    ',
			      -relief => 'raised',
			      -foreground => $ACTIVE_COLOR,
			      -command => sub{&$switch_func($eval_index)}
			      )->grid($options_buttons[$graph_index],
				      -sticky => 'nsew');
    my $window_pane = $window->Scrolled("Pane", 
					-scrollbars => 'e',
					-sticky => 'nswe',
					);
    my $display_frame = $window_pane->Frame(-relief => 'ridge',
					    -borderwidth => 2)->pack(-fill => 'both',
								     -expand => 1);
    $window_pane->pack(-expand => 1, -fill => 'both');
    $options_frames[$eval_index] = $display_frame->Frame()->pack();
    $options_frames[$graph_index] = $display_frame->Frame();
    #Fill eval options frame here
    foreach my $level (sort {$a cmp $b} keys %Display){
	my $level_display  = $level;
	$level_display =~ s/_/ /g;
	my $frame = $options_frames[$eval_index]->Frame(-relief => 'raised',
						       -borderwidth => 1
						       )->grid("-",-sticky => 'ew');
	$frame->Label(-text=> "$t$level_display Level",
		      -justify => 'left')->pack(-side => "left",
						-anchor => 'w');
	my $opt_lev = $Display{$level};
	foreach my $catagory (sort {$a cmp $b} keys %$opt_lev){	    
	    my $catagory_display  = $catagory;
	    $catagory_display =~ s/_/ /g;
	    my $frame = $options_frames[$eval_index]->Frame(-relief => 'ridge',
							    -borderwidth => 1
							    )->grid("-",-sticky => 'ew');
	    $frame->Label(-text=> "$t$t$catagory_display",
			  -justify => 'left')->pack(-side => "left",
						    -anchor => 'w');
	    
	    my $opt_cat = $$opt_lev{$catagory};
	    foreach my $option (sort {$a cmp $b} keys %$opt_cat){		
		my $option_display  = $option;
		$option_display =~ s/_/ /g;
		$options_frames[$eval_index]->Label(-text => "$t$t$t$option_display")
		    ->grid($options_frames[$eval_index]->
			   Checkbutton(-text => "",
				       -variable => 
				       \$Display{$level}{$catagory}{$option}),
			   -sticky => 'nsw');
	    }
	}
    }
    #Fill graph options frame here;
    my $res_frame = $options_frames[$graph_index]->Frame()->pack(-side => 'top',
								 -fill => 'both',
								 -expand => 1);
    $res_frame->Label(-text => "Select Split Type:    ")->pack(-side => 'top');
    my @split_types = Eval::get_graph_x_types();
    my $splitlb = $res_frame->Listbox(-height => $#split_types + 1)->pack();
    foreach my $split (@split_types){
	$splitlb->insert('end',$split);
    }
    my $split_type = $splitlb->get(0);
    $splitlb->selectionSet(0);
    $res_frame->Label(-text => "Select Resolution:    ")->pack(-side => 'top');
    my $uniform_type = 0;
    my $user_defined_type = 1;
    my $graph_split_type = $uniform_type;
    my $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					       -fill => 'x');
    $temp_frame->Radiobutton(-text => "Uniform Resolution",
			     -value => $uniform_type,
			     -variable => \$graph_split_type)->pack(-side => 'left');
    my $get_num_bars = 1;
    my $get_min_x = 0;
    my $get_bar_width = 0;
    my $num_bars = 10;
    my $min_x;
    my $bar_width;
    my $buffer = "    ";
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    my $num_bars_check = $temp_frame->Checkbutton(-text => "Number of Bins",
						  -variable => \$get_num_bars
						  )->pack(-side => 'left');
    my $num_bars_entry = $temp_frame->Entry(-state => 'normal',
					    -textvariable => \$num_bars
					    )->pack(-side => 'right');
    $num_bars_check->bind("<Button-1>",sub {
	if($get_num_bars){
	    $num_bars_entry->configure(state => 'normal',foreground => $ACTIVE_COLOR);
	}
	else{
	    $num_bars_entry->configure(state => 'disabled',foreground => $INACTIVE_COLOR);
	}
    });
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    my $min_x_check = $temp_frame->Checkbutton(-text => "Minimum X Value",
					       -variable => \$get_min_x
					       )->pack(-side => 'left');
    my $min_x_entry = $temp_frame->Entry(-state => 'disabled',
					 -textvariable => \$min_x
					 )->pack(-side => 'right');
    $min_x_check->bind("<Button-1>",sub {
	if($get_min_x){
	    $min_x_entry->configure(state => 'normal',foreground => $ACTIVE_COLOR);
	}
	else{
	    $min_x_entry->configure(state => 'disabled',foreground => $INACTIVE_COLOR);
	}
    });
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    my $bar_width_check = $temp_frame->Checkbutton(-text => "Bin Width",
						   -variable => \$get_bar_width
						   )->pack(-side => 'left');
    my $bar_width_entry = $temp_frame->Entry(-state => 'disabled',
					     -textvariable => \$bar_width
					     )->pack(-side => 'right');
    $bar_width_check->bind("<Button-1>",sub {
	if($get_bar_width){
	    $bar_width_entry->configure(state => 'normal',foreground => $ACTIVE_COLOR);
	}
	else{
	    $bar_width_entry->configure(state => 'disabled',foreground => $INACTIVE_COLOR);
	}
    });
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');    
    $temp_frame->Radiobutton(-text => "User Defined Resolution",
			     -value => $user_defined_type,
			     -variable => \$graph_split_type)->pack(-side => 'left');
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');    
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    $temp_frame->Label(-text => "From")->pack(-side => 'left');
    my $from;
    $temp_frame->Entry(-textvariable => \$from)->pack(-side => 'right');
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    my $to;
    $temp_frame->Label(-text => "To")->pack(-side => 'left');
    $temp_frame->Entry(-textvariable => \$to)->pack(-side => 'right');
    $temp_frame = $res_frame->Frame()->pack(-side => 'top',
					    -fill => 'x');    
    $temp_frame->Label(-text => $buffer)->pack(-side => 'left');
    my $reslb = $res_frame->Listbox();
    $temp_frame->Button(-text => 'Add',
			command => sub{
			    unless(($from =~ /^\d+\.\d+$/) ||
				   ($from =~ /^\d+$/) ||
				   ($from =~ /^\.\d+$/)){
				error_func("\"From\" value is not a number.  This range".
					   " will not be added.");
				return;
			    }
			    unless(($to =~ /^\d+\.\d+$/) ||
				   ($to =~ /^\.\d+$/) ||
				   ($to =~ /^\d+$/)){
				error_func("\"To\" value is not a number.  This range".
					   " will not be added.");
				return;
			    }
			    unless($to > $from){
				error_func("\"To\" value must be greater than \"From\"".
					   "value.  This range will not be added.");
			    }
			    my @all = $reslb->get(0,'end');
			    if($#all == -1){
				$reslb->insert('end',"$from -> $to");
				$from = $to;
				$to = "";
				return;
			    }
			    my @ranges;
			    foreach my $range (@all){
				$range =~ /(\S+) -> (\S+)/;
				push @ranges, {start => $1,
					       stop => $2};
			    }
			    $reslb->delete(0,'end');
			    push @ranges, {start => $from,
					   stop => $to};
			    @ranges = sort {$$a{start} <=> $$b{start}} @ranges;
			    my $last_end = -1;
			    my $error = 0;
			    foreach my $range (@ranges){
				if($$range{start} < $last_end){
				    unless($error){
					error_func("New range overlaps existing ranges.".
						   "  The ranges will be adjusted.");
					$error = 1;
				    }
				    $$range{start} = $last_end;
				}
				$last_end = $$range{stop};
				$reslb->insert('end',"$$range{start} -> $$range{stop}"); 
			    }
			    $from = $to;
			    $to = "";
			})->pack(-side => 'left');
    $temp_frame->Button(-text => 'Remove',
			-command => sub{
			    my @sel = $reslb->curselection();
			    while($#sel >= 0){
				$reslb->delete($sel[0]);
				@sel = $reslb->curselection();
			    }
			})->pack(-side => 'right');
    $reslb->pack(-fill => 'x');
    my $set_graph_res_sub = sub {
	my ($graph_type) = @_;
	if($graph_split_type == $uniform_type){
	    if(defined($Graph_Resolution{$graph_type}{user})){
		delete($Graph_Resolution{$graph_type}{user});
	    }
	    $Graph_Resolution{$graph_type}{uniform} = {count => '',
						       min => '',
						       size => ''};
	    if($get_num_bars){
		$Graph_Resolution{$graph_type}{uniform}{count} = $num_bars;
	    }
	    if($get_min_x){
		$Graph_Resolution{$graph_type}{uniform}{min} = $min_x;
	    }
	    if($get_bar_width){
		$Graph_Resolution{$graph_type}{uniform}{size} = $bar_width;
	    }
	}
	elsif($graph_split_type == $user_defined_type){
	    if(defined($Graph_Resolution{$graph_type}{uniform})){
		delete($Graph_Resolution{$graph_type}{uniform});
	    }
	    my @ranges;
	    $Graph_Resolution{$graph_type}{user} = \@ranges;
	    my @temp = $reslb->get(0,'end');	    
	    my $last_stop = $temp[0];
	    $last_stop =~ /(\S+) -> (\S+)/;
	    $last_stop = $1;
	    my $bad = 0;
	    $reslb->delete(0,'end');
	    foreach my $range (@temp){
		$range =~ /(\S+) -> (\S+)/;
		my $min = $1;
		my $max = $2;
		unless($min == $last_stop){
		    push @ranges, {min => $last_stop,
				   max => $min};
		    $reslb->insert('end',"$last_stop -> $min");
		    $bad = 1;
		}
		$last_stop = $max;
		push @ranges, {min => $min,
			       max => $max};
		$reslb->insert('end',"$min -> $max");
	    }
	    if($bad){
		error_func("No gaps allowed between bins.  New bins have been added to".
			   " fill the gaps.");
	    }
	}
	else{
	    #EVAN error
	}	
    };
    my $load_graph_res_sub = sub {
	my ($graph_type) = @_;
	if(defined($Graph_Resolution{$graph_type}{uniform})){
	    $graph_split_type = $uniform_type;
	    $num_bars = $Graph_Resolution{$graph_type}{uniform}{count};
	    if($num_bars eq ''){
		$get_num_bars = 0;
		$num_bars_entry->configure(state => 'disabled',
					   foreground => $INACTIVE_COLOR);
	    }
	    else{
		$get_num_bars = 1;
		$num_bars_entry->configure(state => 'normal',
					   foreground => $ACTIVE_COLOR);
	    }
	    $min_x = $Graph_Resolution{$graph_type}{uniform}{min};
	    if($min_x eq ''){
		$get_min_x = 0;
		$min_x_entry->configure(state => 'disabled',
					   foreground => $INACTIVE_COLOR);
	    }
	    else{
		$get_min_x = 1;
		$min_x_entry->configure(state => 'normal',
					foreground => $ACTIVE_COLOR);
	    }
	    $bar_width = $Graph_Resolution{$graph_type}{uniform}{size};
	    if($bar_width eq ''){
		$get_bar_width = 0;
		$bar_width_entry->configure(state => 'disabled',
					    foreground => $INACTIVE_COLOR);
	    }
	    else{
		$get_bar_width = 1;
		$bar_width_entry->configure(state => 'normal',
					    foreground => $ACTIVE_COLOR);
	    }
	    $reslb->delete(0,'end');
	}
	elsif(defined($Graph_Resolution{$graph_type}{user})){
	    $graph_split_type = $user_defined_type;
	    my $ranges = $Graph_Resolution{$graph_type}{user};
	    $reslb->delete(0,'end');
	    foreach my $range (@$ranges){
		$reslb->insert('end',"$$range{min} -> $$range{max}");
	    }
	    $num_bars = '';
	    $min_x = '';
	    $bar_width = '';
	    $get_num_bars = 0;
	    $get_min_x = 0;
	    $get_bar_width = 0;
	    $num_bars_entry->configure(state => 'disabled',
				       foreground => $INACTIVE_COLOR);
	    $min_x_entry->configure(state => 'disabled',
				    foreground => $INACTIVE_COLOR);
	    $bar_width_entry->configure(state => 'disabled',
					foreground => $INACTIVE_COLOR);
	}
	else{
	    #EVAN error
	}	
    };
    &$load_graph_res_sub($split_type);
    $splitlb->bind("<Button-1>",sub {&$set_graph_res_sub($split_type);
				     my @sel = $splitlb->curselection();
				     $split_type = $splitlb->get($sel[0]);
				     &$load_graph_res_sub($split_type)});
    my %default = get_default_graph_resolution();
    $button_frame->Button(-text => "Save Changes",
			  -command => sub{
			      &$set_graph_res_sub($split_type);
			      save_options_file()})->pack(-side => 'left');
    $button_frame->Button(-text => "Close",
			  -command => sub {$window->destroy})->pack(-side => 'right');   
}

__END__
