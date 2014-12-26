#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my ($samples,$data_type,$suffix1,$suffix2,$help);
GetOptions('f=s' => \$samples,
	   'd=s' => \$data_type,
	   'a=s' => \$suffix1,
	   'b=s' => \$suffix2,
	   'h' => \$help);

my $oyts=$ENV{'LAOTS'};

my @samples = ();
if ( -e $samples ) {
    open my $in, '<', $samples;
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);
	push(@samples, $f[0]);
    }
} else {
    push(@samples, $samples);
}

print qq(sample\tdat_type\tuniq\tdups\tdup_percent\n);
foreach my $s (@samples) {
    my $reads1 = "$oyts/$s/$data_type/${s}_${data_type}.$suffix1.R1.fq";
    my $reads2 = "$oyts/$s/$data_type/${s}_${data_type}.$suffix1.$suffix2.R1.fq";
    if ( -e $reads1 ) {
	my $wc1 = `wc -l $reads1 | cut -f 1 -d ' '`;
	my $wc2 = `wc -l $reads2 | cut -f 1 -d ' '`;
	chomp($wc1); chomp($wc2);
	$wc1 = $wc1/4;
	$wc2 = $wc2/4;
	my $dups = $wc1 - $wc2;
	my $perc = sprintf("%.2f", 100*($dups/$wc1));
	print qq($s\t$data_type\t$wc2\t$dups\t$perc\n);
    }    
}

exit;
