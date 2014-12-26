#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Getopt::Long;

my ($file,$ref,$suffix,$work_dir,$log,$help);
my $dat_type = 'MG';
my $win_lng = 50;
GetOptions('f|file=s' => \$file,
	   'r|ref=s' => \$ref,
	   'd|data_type=s' => \$dat_type,
	   'a|work_dir=s' => \$work_dir,
	   'w|window=s' => \$win_lng,
	   'z|suffix=s' => \$suffix,
	   'l|log' => \$log);

die if (!$file || !$ref);
&loadEnv('/mnt/nfs/projects/ecosystem_biology/esb.shortcuts');
&usage if ($help || !$file || !$ref);

my $samples = &parseSamples($file);

foreach my $s (sort keys %{$samples}) {
    print qq(Processing $s.\n);
    my $sequences = &loadSeqs($ref); # this is inefficient.
                                     # must add way of purging the cov data on each pass

    my $cov = "$ENV{'LAOTS'}/$s/$dat_type/$work_dir/$s\_$dat_type.$suffix.cov";
    my $ofile = $cov;
    if ($log) {
	$ofile =~ s/\.cov$/\.win${win_lng}bp\.log10\.cov/;
    } else {
	$ofile =~ s/\.cov$/\.win${win_lng}bp\.cov/;
    }

    open my $in, '<', $cov;
    print qq(Gathering coverage data for $s.\n);
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);    
	$$sequences{$f[0]}{'position'}{$f[1]} = $f[2];
    }
    close($in);

    print qq(Starting coverage calculations for $s.\n);
    open my $out, '>', $ofile;
    foreach my $seq (sort {$$sequences{$b}{'length'} <=> $$sequences{$a}{'length'}} keys %{$sequences}) {
	my $start = 1;
	my $end = $win_lng;
	my $cum = 0;
	my $ave;
    
	foreach my $pos (sort {$a <=> $b} keys %{$$sequences{$seq}{'position'}}) {
	    if ($pos > $end) {
		$ave = sprintf("%.0f",$cum/$win_lng);
		if ($log) {
		    $ave = &log10($ave + 1);
		    $ave = sprintf("%.4f", $ave);
		} else {
		    $ave = sprintf("%.2f", $ave);
		}
		print $out qq($seq $start $end $ave\n);
		
		# reset values
		$start = $start + $win_lng;
		$end = $start + $win_lng - 1;
		$cum = $$sequences{$seq}{'position'}{$pos};
		
	    } else {
		$cum += $$sequences{$seq}{'position'}{$pos};
		
		if ($pos == $$sequences{$seq}{'length'}) {
		    # print the last window for each sequence
		    $ave = sprintf("%.0f",$cum/($pos-$start+1));
		    if ($log) {
			$ave = &log10($ave + 1);
			$ave = sprintf("%.4f", $ave);
		    } else {
			$ave = sprintf("%.2f", $ave);
		    }
		    print $out qq($seq $start $pos $ave\n);	
		}
	    }
	}
    }
}

exit;
    
###############
# subroutines #
###############
sub loadEnv {
    my $f = shift;
    my $perl_command = "perl -MData::Dumper -e 'print Dumper(\\\%ENV)';";
    my $source_line = ". $f 1>&2";
    my %tmp = %{eval('my ' . `$source_line\n$perl_command`)};
    $ENV{$_} = $tmp{$_} for (keys %tmp);
    return 0;
}

sub loadSeqs {
    my $f = shift;
    my %hash = ();
    my $iseq = Bio::SeqIO->new(-file => $f,
			       -format => 'fasta');
    while (my $seqobj = $iseq->next_seq()) {
	$hash{$seqobj->display_id()}{'length'} = $seqobj->length();
	for (my $i=1;$i<=$seqobj->length();$i++) {
	    $hash{$seqobj->display_id()}{'position'}{$i} = 0;
	}
    }
    return \%hash;
}

sub parseSamples {
    my $f = shift;
    my %hash = ();
    if (-e $f) {
	open my $in, '<', $f;
	while (<$in>) {
	    next if (/^\#/);
	    chomp;
	    my @s = split(/\t/);
	    $hash{$s[0]}{'date'} = $s[1];
	}
    } else { $hash{$f} = 1; } # this will give an error, since storing the matrix relies on dates. fix above
    return \%hash;
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
