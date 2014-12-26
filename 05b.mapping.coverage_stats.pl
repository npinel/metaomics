#!/usr/bin/perl -w

=head1 NAME

  05b.mapping.coverage_stats.pl - one line description

=head1 SYNOPSIS

    Brief summary of what this script does

=head1 DESCRIPTION

    Describe what is the purpose and characteristics of this script

=head1 USAGE

=head1 AUTHOR - Nicolas Pinel (March 19, 2013)

=head1 PENDING

    * NOTE: for some unidentified reason, this script is not working; it runs silently,
      but the expected output file is not produced
    * To implement the option to speficy bam files directly rather than assume their path

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Statistics::Descriptive;

$Term::ANSIColor::AUTORESET = 1; # to reset color to defaul on the next print statement 

my ($ref,$ref_desc,$file,$mapper,$work_dir,$suffix,$verbose,$help);
my $min_cov = 1;
my $min_mapq = 5;
my $dat_type = 'MG';
GetOptions('r|ref=s'      => \$ref,
	   'i|ref_desc=s' => \$ref_desc,
 	   'd|dat_type=s' => \$dat_type,
	   'f|file=s'     => \$file,
	   'm|mapper=s'   => \$mapper,
	   'w|dir=s'      => \$work_dir,
	   'c|min_cov=s'  => \$min_cov,
	   'q|min_mapq=s' => \$min_mapq,
	   'z|suffix=s'   => \$suffix,
	   'v|verbose'    => \$verbose,
	   'h|help'       => \$help);

&loadEnv('/mnt/nfs/projects/ecosystem_biology/esb.shortcuts');
&usage if ($help || !$file || !$ref);

# system("module load SAMtools/0.1.18-ictce-5.3.0");
my $samtools = "$ENV{'TOOLS'}/samtools-0.1.19/bin/samtools";
my $samples = &parseSamples($file);

my ($contig_lng,$tot_ref_lng) = &getRefStats($ref);
    # $contig_lng is a hash_ref with contig_id => contig_length
    # $tot_ref_lng is a scalar with the value of the total seq length

my %stats = ();
foreach my $s (sort keys %{$samples}) {
    print qq(processing $s\n) if ($verbose);
    my $bam = "$ENV{'LAOTS'}/$s/$dat_type/$work_dir/$s\_$dat_type.$mapper.${ref_desc}_x_$suffix.bam";

    unless (-e $bam) {
	print BOLD RED 'ERROR:', RESET, qq(\tThe bam file $bam could not be found. Moving to the next sample...\n);
    }
    
    my $glb_cov_tot = 0;
    my $glb_cov = 0;
    my $glb_cov_lng = 0;
    my %cov = ();
    my $cov  = $bam;
    $cov =~ s/\.bam$/\.cov/;
    if ( -e $cov) {
	print qq(found coverage file $cov. proceeding with tally...\n);
    } else {
	print qq(could not find coverage file. will compute for $bam.\n);
#	my $cmd = "samtools mpileup -q $min_mapq -f $ref $bam 2>/dev/null | cut -f 1,2,4 > $cov";
	my $cmd = "$samtools depth -Q $min_mapq $bam 2>/dev/null > $cov";
	system($cmd);
    }

     open my $in, '<', "$cov";
    my $new_contig = '';
    my $old_contig = '';
     while (<$in>) {
 	chomp;
 	my @set = split(/\t/);
	$new_contig = $set[0];
 	## gather global stats
# 	$glb_cov_tot += $set[2];
# 	$cov{$set[0]}{'cov_tot'} = $cov{$set[0]}{'cov_tot'} ? $cov{$set[0]}{'cov_tot'} + $set[2] : $set[2];
# 	$stats{$set[0]}{$$samples{$s}{'date'}}{'cov_tot'} = $stats{$set[0]}{$$samples{$s}{'date'}}{'cov_tot'} ? $stats{$set[0]}{$$samples{$s}{'date'}}{'cov_tot'} + $set[2] : $set[2];
 	## gather stats above cov cutoff
 	if ($set[2] >= $min_cov) {
#	    $glb_cov += $set[2];
# 	    $glb_cov_lng++;
# 	    $cov{$set[0]}{'cov'} = $cov{$set[0]}{'cov'} ? $cov{$set[0]}{'cov'} + $set[2] : $set[2];
# 	    $cov{$set[0]}{'lng'} = $cov{$set[0]}{'lng'} ? $cov{$set[0]}{'lng'} + 1 : 1;
 	    $stats{$set[0]}{$$samples{$s}{'date'}}{'cov'} = $stats{$set[0]}{$$samples{$s}{'date'}}{'cov'} ? $stats{$set[0]}{$$samples{$s}{'date'}}{'cov'} + $set[2] : $set[2];
 	    $stats{$set[0]}{$$samples{$s}{'date'}}{'lng'} = $stats{$set[0]}{$$samples{$s}{'date'}}{'lng'} ? $stats{$set[0]}{$$samples{$s}{'date'}}{'lng'} + 1 : 1;
 	}
	print qq(\t$new_contig\n) if ($verbose && ($old_contig ne $new_contig));
	$old_contig = $new_contig;
     }
     ## output stats for each file
#    if ($verbose) {
# 	my $glb_cov_tot_ave = sprintf("%.2f", ($glb_cov_tot/$tot_ref_lng));
# 	my $glb_cov_ave = sprintf("%.2f", ($glb_cov/$glb_cov_lng));
# 	my $glb_cov_frc = sprintf("%.2f", ($glb_cov_lng/$tot_ref_lng));

# 	print qq(\# Stats\: $s\nglobal mean coverage: $glb_cov_tot_ave\nabove $min_cov x\n\tmean coverage\: $glb_cov_ave\n\tcovered fraction\: $glb_cov_frc\n);
# 	print qq(\# Coverage by contig:\n);
# 	foreach my $contig (sort keys %cov) {
# 	    my $cov_tot_ave = sprintf("%.2f", ($cov{$contig}{'cov_tot'}/$$contig_lng{$contig}));
# 	    my $cov_ave = ($cov{$contig}{'lng'} && $cov{$contig}{'lng'} > 0) ? sprintf("%.2f", ($cov{$contig}{'cov'}/$cov{$contig}{'lng'})) : 0;
# 	    my $cov_frc = ($cov{$contig}{'lng'} && $cov{$contig}{'lng'} > 0) ? sprintf("%.2f", ($cov{$contig}{'lng'}/$$contig_lng{$contig})) : 0;
# 	    print qq(\t\* contig\:$contig\n\t\tmean coverage: $cov_tot_ave\n),
# 	    qq(\t\tabove $min_cov x\n\t\t\tmean coverage\: $cov_ave\n\t\t\tcovered fraction\: $cov_frc\n);
# 	}
#     }
}

chomp(my $date = `date +%Y%m%d`);
open my $out, '>', "$ENV{'LAOTS'}/coverage_matrix.$date.$ref_desc.txt";

my $str = "date\t".join("\t", sort keys %stats);
print $out $str,"\n";
foreach my $s (sort {$$samples{$a}{'date'} cmp $$samples{$b}{'date'}} keys %{$samples}) {
    $str = $$samples{$s}{'date'};
    foreach my $ctg (sort keys %stats) {
	my $ave = ($stats{$ctg}{$$samples{$s}{'date'}}{'lng'} && $stats{$ctg}{$$samples{$s}{'date'}}{'lng'} > 0) ? sprintf("%.2f", ($stats{$ctg}{$$samples{$s}{'date'}}{'cov'}/$stats{$ctg}{$$samples{$s}{'date'}}{'lng'})) : 0;
	$str .= "\t$ave";
#	my $cov_frc = ($cov{$contig}{'lng'} && $cov{$contig}{'lng'} > 0) ? sprintf("%.2f", ($cov{$contig}{'lng'}/$$contig_lng{$contig})) : 0;
    }
    print $out $str,"\n";
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
    } else { $hash{$f} = 1; } # this will give an error, since storing the matrix relies on dates. fix above

    return \%hash;
}

sub getRefStats {
    my $f = shift;

    my %hash = ();
    my $total = 0;
    my $seqin = Bio::SeqIO->new('-file' => $f,
				'-format' => 'fasta');

    while (my $seqobj = $seqin->next_seq()) {
	$hash{$seqobj->display_id()} = $seqobj->length();
	$total += $seqobj->length();
    }

    return \%hash, $total;
}
