#! /usr/bin/env perl

$FIRST_READ_P = 0x0040;  # 64
$SECOND_READ_P = 0x0080; # 128
$QUERY_MINUSP = 0x0010;  # 16


while (defined($line = <>)) {
    @fields = split /\t/,$line;

    $acc = $fields[0];

    $flags = $fields[1];
    if ($flags & $FIRST_READ_P) {
	$acc .= "/1";
    } elsif ($flags & $SECOND_READ_P) {
	$acc .= "/2";
    }
    if ($flags & $QUERY_MINUSP) {
	$plusp = 0;
    } else {
	$plusp = 1;
    }

    $hiti = get_aux("HI:i",\@fields);
    $acc .= "." . $hiti;

    $nm_score = get_aux("NM:i",\@fields);
    $xs_score = get_aux("XS:A",\@fields);

    $chr = $fields[2];
    $chrpos = $fields[3];	# Excludes soft clip
    $cigar = $fields[5];
    $read = $fields[9];
    $quality_string = $fields[10];

    process_exons($acc,$cigar,$plusp,$read,$quality_string,$chr,$chrpos,
 	          $nm_score,$xs_score);
}

exit;


sub print_exon {
    my ($acc, $plusp, $read, $quality_string,
	$readpos, $npos, $chr, $chrpos, $nm_score, $xs_score) = @_;

    print ">$acc" . "@";

    if ($plusp == 1) {
	printf "%d",$readpos + 1;
    } else {
	printf "-%d",$readpos + 1;
    }

    printf " $chr:%u..%u M",$chrpos,$chrpos + $npos - 1;
    if (defined($nm_score)) {
	printf " NM:i:%d",$nm_score;
    }
    if (defined($xs_score)) {
	printf " XS:A:%d",$xs_score;
    }
    print "\n";

    print substr($read,$readpos,$npos) . "\n";

    if ($quality_string ne "*") {
	print substr($quality_string,$readpos,$npos) . "\n";
    }

    return;
}

sub print_deletion {
    return;
}

sub print_insertion {
    return;
}

sub print_low_softclip {
    return;
}

sub print_high_softclip {
    return;
}


sub process_exons {
    my ($acc, $cigar, $plusp, $read, $quality_string, $chr, $chrpos,
	$nm_score, $xs_score) = @_;
    my $readpos = 0;
    
    ($npos,$type,$rest) = $cigar =~ /^(\d+)(\D)(.*)/;

    if (!defined($type)) {
	die "Cannot parse cigar $cigar";

    } elsif ($type eq "M") {
	print_exon($acc,$plusp,$read,$quality_string,
		   $readpos,$npos,$chr,$chrpos,$nm_score,$xs_score);
	$readpos += $npos;
	$chrpos += $npos;

    } elsif ($type eq "N") {
	$chrpos += $npos;
	
    } elsif ($type eq "D") {
	print_deletion($part++);
	$chrpos += $npos;
	
    } elsif ($type eq "I") {
	print_insertion($part++);
	$readpos += $npos;
	
    } elsif ($type eq "S") {
	print_low_softclip($part++);
	$readpos += $npos;
    }
    
    $cigar = $rest;

    while ($cigar =~ /\S/) {
	($npos,$type,$rest) = $cigar =~ /^(\d+)(\D)(.*)/;
	if (!defined($type)) {
	    die "Cannot parse cigar $cigar";

	} elsif ($type eq "M") {
	    print_exon($acc,$plusp,$read,$quality_string,
		       $readpos,$npos,$chr,$chrpos,$nm_score,$xs_score);
	    $readpos += $npos;
	    $chrpos += $npos;

	} elsif ($type eq "N") {
	    $chrpos += $npos;

	} elsif ($type eq "D") {
	    print_deletion($part++);
	    $chrpos += $npos;

	} elsif ($type eq "I") {
	    print_insertion($part++);
	    $readpos += $npos;

	} elsif ($type eq "S") {
	    print_high_softclip($part++);
	    $readpos += $npos;
	}

	$cigar = $rest;
    }
    
    return;
}


sub get_aux {
    my ($desired_field, $fields) = @_;
    my $fieldi = 11;

    while ($fieldi <= $#{$fields}) {
	if ($ {$fields}[$fieldi] =~ /^$desired_field:(\S+)/) {
	    return $1;
	}
	$fieldi++;
    }
    return;
}


