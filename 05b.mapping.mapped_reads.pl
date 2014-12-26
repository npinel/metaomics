#!/usr/bin/perl -w

=head1 NAME

  05b.mapping.mapped_reads.pl - one line description

=head1 SYNOPSIS

    Brief summary of what this script does

=head1 DESCRIPTION

    Describe what is the purpose and characteristics of this script

=head1 USAGE

=head1 AUTHOR - Nicolas Pinel (1 May 2013)

=head1 PENDING

    * To implement the option to speficy bam files directly rather than assume their path

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Statistics::Descriptive;

$Term::ANSIColor::AUTORESET = 1; # to reset color to defaul on the next print statement 

my ($ref_desc,$file,$bam,$verbose,$help);
my $min_mapq = 5;
my $dat_type = 'MG';
GetOptions('i|ref_desc=s' => \$ref_desc,
	   'd|dat_type=s' => \$dat_type,
	   'b|bam=s'      => \$bam,
	   'f|file=s'     => \$file,
	   'q|min_mapq=s' => \$min_mapq,
	   'v|verbose'    => \$verbose,
	   'h|help'       => \$help);

&loadEnv('/mnt/nfs/projects/ecosystem_biology/esb.shortcuts');
&usage if ($help || !$file);

my $samples = &parseSamples($file);

foreach my $s (@{$samples}) {
    $bam = "$ENV{'LAOTS'}/$s/$dat_type/$s\_$dat_type.mapped_$ref_desc.bam";
    unless (-e $bam) {
	print BOLD RED 'ERROR:', RESET, qq(\tThe bam file $bam could not be found. Moving to the next sample...\n);
    }

    chomp(my $mapped_tot = `samtools view $bam | wc -l`);
    chomp(my $mapped_hmq = `samtools view -q 5 $bam | wc -l`);
    print qq($s\t$mapped_tot\t$mapped_hmq\n);
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

    my @array = ();
    if (-e $f) {
	open my $in, '<', $f;
	while (<$in>) {
	    next if (/^\#/);
	    my @s = split(/\t/);
	    push(@array, $s[0]);
	}
    } else { @array = $f; }

    return \@array;
}
