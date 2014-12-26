#!/usr/bin/perl -w

=head1 NAME

=head1 SYNOPSIS

    Brief summary of what this script does

=head1 DESCRIPTION

    Describe what is the purpose and characteristics of this script

=head1 USAGE

=head1 AUTHOR - Nicolas Pinel (20 September 2013)

=cut

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

$Term::ANSIColor::AUTORESET = 1; # to reset color to defaul on the next print statement 

my ($sample_list,$data_type,$targets,$annotation,$include_rna,$dat_dir_mg,$dat_dir_mt,$file_info_mg,$file_info_mt,$log,$pseudo_zero,$verbose,$help);
my $out_fmt = 'cmonkey';
my $metric = 'cov';
my $na_str = 'NA';
my $pzv = 0.00001; # pseudo-zero value
GetOptions('a|annotation=s' => \$annotation,
	   'o|out_format=s' => \$out_fmt,
	   't|type=s'       => \$data_type,
	   'r|include_rna'  => \$include_rna,
	   's|samples=s'    => \$sample_list,
	   'd|dir=s'        => \$dat_dir_mg,
	   'x=s'            => \$dat_dir_mt,
	   'y=s'            => \$file_info_mt,
	   'i|info=s'       => \$file_info_mg,
	   'm|metric=s'     => \$metric,
	   'n|na_str=s'     => \$na_str,
	   'l|log'          => \$log,
	   'z|pseudo_zero'  => \$pseudo_zero,
	   'v|verbose'      => \$verbose,
	   'h|help'         => \$help);

&loadEnv('/mnt/nfs/projects/ecosystem_biology/esb.shortcuts');
&usage if ($help || !$sample_list || !$data_type);
$file_info_mg =~ s/^\.//;
$file_info_mg =~ s/\.$//;
$file_info_mt =~ s/^\.//;
$file_info_mt =~ s/\.$//;

my $samples = &parseSamples($sample_list);

chomp(my @features = `cut -f 9 $annotation`);
my %features = ();

foreach my $s (keys %{$samples}) {
    my %data = ('MG' => "$ENV{'LAOTS'}/$s/MG/$dat_dir_mg/${s}.$file_info_mg.isoforms.fpkm_tracking",
		'MT' => "$ENV{'LAOTS'}/$s/MT/$dat_dir_mt/${s}.$file_info_mt.isoforms.fpkm_tracking");

    if ($data_type =~ /R/) {
	unless ( (-e $data{'MG'}) && (-e $data{'MT'}) ) {
	    warn("cannot find $s\n");
	    delete($$samples{$s});
	}
	foreach my $data (keys %data) {
	    open my $in, '<', $data{$data};
	    while (<$in>) {
		next if (/^tracking_id/);
		chomp;
		my @fields = split(/\t/);
		next if ($fields[0] =~ /rna/ && !$include_rna);		
		$features{$fields[0]}{$$samples{$s}{'date'}}{$data} = $metric =~ /cov/i ? $fields[8] : $fields[9]; # to print in cMonkey format
	    }
	}
    } else {
	    open my $in, '<', $data{$data_type};
	    while (<$in>) {
		next if (/^tracking_id/);
		chomp;
		my @fields = split(/\t/);
		next if ($fields[0] =~ /rna/ && !$include_rna);		
		$features{$fields[0]}{$$samples{$s}{'date'}}{$data_type} = $metric =~ /cov/i ? $fields[8] : $fields[9]; # to print in cMonkey format
	    }	
    }

}

my $str;
if ($out_fmt =~ /cmonkey/i) {
    $str = "feat_x_date";
    foreach my $k (sort {$$samples{$a}{'date'} cmp $$samples{$b}{'date'}} keys %{$samples}) {
	$str .= "\t$$samples{$k}{'date'}";
    }
    print $str,"\n";

    foreach my $f (sort keys %features) {
	$str = $f;
	foreach my $k (sort {$$samples{$a}{'date'} cmp $$samples{$b}{'date'}} keys %{$samples}) {
	    if ($data_type =~ /RT/i) {
		my $mg = $features{$f}{$$samples{$k}{'date'}}{'MG'};
		my $mt = $features{$f}{$$samples{$k}{'date'}}{'MT'};
		my $ratio = $mg == 0 ? $na_str : ($mt == 0 ? ($pseudo_zero ? $pzv : $na_str) : 
						  ( $data_type =~ /R/ ? ($log ? sprintf("%.4f", &log2($mt/$mg)) : sprintf("%.4f", $mt/$mg)) :
						    sprintf("%.4f", $features{$f}{$$samples{$k}{'date'}}{$data_type})));
		$str .= "\t$ratio";
	    } else {
		$str .= "\t$features{$f}{$$samples{$k}{'date'}}{$data_type}";
	    }
	}
	print $str,"\n";
    }
    
} else {
	
    $str = "date_x_feat\t".join("\t", sort keys %features);
    print $str,"\n";	
    foreach my $k (sort {$$samples{$a}{'date'} cmp $$samples{$b}{'date'}} keys %{$samples}) {
	$str = $$samples{$k}{'date'};
	foreach my $f (sort keys %features) {
	    if ($data_type =~ /RT/i) {
		my $mg = $features{$f}{$$samples{$k}{'date'}}{'MG'};
		my $mt = $features{$f}{$$samples{$k}{'date'}}{'MT'};
		my $ratio = $mg == 0 ? $na_str : ($mt == 0 ? ($pseudo_zero ? $pzv : $na_str) :
						  ( $data_type =~ /R/ ? ($log ? sprintf("%.4f", &log2($mt/$mg)) : sprintf("%.4f", $mt/$mg)) :
						    sprintf("%.4f", $features{$f}{$$samples{$k}{'date'}}{$data_type})));
		$ratio = $na_str if ($ratio == 0);
		$str .= "\t$ratio";
	    } else {
		$str .= "\t$features{$f}{$$samples{$k}{'date'}}{$data_type}";
	    }
	}
	print $str,"\n";
    }
}
exit;

###############
# subroutines #
###############
sub usage {
    system("perldoc $0");
    exit 1;
}

sub loadEnv {
    my $f = shift;
    my $perl_command = "perl -MData::Dumper -e 'print Dumper(\\\%ENV)';";
    my $source_line = ". $f 1>&2";
    my %tmp = %{eval('my ' . `$source_line\n$perl_command`)};
    $ENV{$_} = $tmp{$_} for (keys %tmp);
    return 0;
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
    } else { $hash{$f} = 1; }

    return \%hash;
}

sub log2 {
    my $n = shift;
    return log($n)/log(2);
}
