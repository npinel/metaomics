#!/usr/bin/perl -w

=head1 NAME

    03c.profile.extract_reads.pl - extract fastq reads from indexed files

=head2 Dependencies

    cdbyank - for extracting the reads from indexed sequence reads files

=head2 Description

    Use this script to extract (subsets of) sequence reads from indexed reads files. Sequence
    read files should be indexed with cdbfasta with default parameter values for keys.

    A list of sequence read keys consists of sequence read ids, but without the sequence 
    delimiters (i.e., '@' in fastq-formatted reads, or '>' in fasta-formatted reads). Keys can
    be obtained, among other ways, from:
        - the output of step 03b (binning by %GC/length)
        - bam files containing mapped reads (filtered with samtools view -F 0x0004)
        - extracting read ids from sequence read files; sequence delimeters must be removed in
          order for the ids to matched the indexing keys

    Currently (25 Apr 2013) this script produces output files. If invoked from within the assembly
    or mapping wrappers, these output files may be deleted in order to minimize storage demands.
    It may also be of interest to integrate this script so as to pipe directly into subsequent steps.

=head2 Input files

    Keys file(s)  - single file if extracting PandaSeq files from R1/R2 files
                  - single file if extracting mapped reads, even if extracting from R1/R2 files
                  - two (paired) files, if extracting reads corresponding to reads from %GC bins,
                    so that both members of the pair get included in to the bin
    Index file(s) - single file if extracting from un-paired reads
                  - two files (*.R1.fq.cidx and *.R2.fq.cidx) whenever paired reads are desired

=head1 USAGE

    ./03c.profile.extract_reads.pl -k1 <keys1> -i1 <reads1> [-i2 <index2>] [-k2 <keys2>] [-t]

        -k1/-keys1    file containing a list of read keys to extract (required)
        -k2/-keys2    second file containing a list of read keys (optional)
        -i1/-index1   cdbfasta-constructed read index file, as input to cdbyank (required)
        -i2/-index2   index file paired to -i1 (optional, but required if paired reads are desired)
        -b/-bam       indicates 'key' file(s) are given as bam alignment files, in which case the keys
                      to use for read extraction will be obtained using samtools view; this option
                      is only available for single key files
        -t/-ps        flag to apply trimming of adaptor indicator; use this option when extracting
                      PandaSeq-joined reads from R1/R2 files, since PandaSeq appends ':NNNNNNNNN' to
                      read identifiers. WARNING: without this flag, keys from PandaSeq files will not
                      be found in regular R1/R2 files
        -e/-exclude   exclude reads from indexes instead of extract
        -c/-cdbyank   exact location of cdbyank, if not on path
        -s/-samtools  exact location of samtools, if not on path
        -v/-verbose   run status messages printed to stderr
        -h/-help      displays this documentation

=head1 AUTHOR - Nicolas Pinel (ISB/LCSB - 25 Apr 2013)

=cut

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

$Term::ANSIColor::AUTORESET = 1; # to reset color to defaul on the next print statement 

my ($keys1,$keys2,$index1,$index2,$trim_adapt,$mode,$bam,$exclude,$help,$verbose);
my $format = 'separate'; # vs inter(leaved)

#system("module load SAMtools/0.1.18-ictce-5.3.0");
chomp(my $cdbyank = `which cdbyank`);
chomp(my $samtools = `which samtools`);
GetOptions('k1|keys1=s'   => \$keys1,
	   'k2|keys2=s'   => \$keys2,
	   'i1|index1=s'  => \$index1,
	   'i2|index2=s'  => \$index2,
	   'f|format=s'   => \$format, ## output for format; separate vs interleaved; NOT YET IMPLEMENTED
	   't|ps'         => \$trim_adapt, ## to indicate whether keys come from PandaSeq
	   'm|mode=s'     => \$mode, ## see documentation
	   'c|cdbyank=s'  => \$cdbyank,
	   'b|bam'        => \$bam,
	   's|samtools=s' => \$samtools,
	   'e|exclude'    => \$exclude,
	   'h|help'       => \$help,
	   'v|verbose'    => \$verbose);

&usage if ($help);
&loadEnv('/mnt/nfs/projects/ecosystem_biology/esb.shortcuts');

$cdbyank = $cdbyank ? $cdbyank : "$ENV{'TOOLS'}/cdbyank";
&checkExecutables();

&checkInputs();

## done with checking, now execute
if ($exclude) {
    &executeExclude();
} else {
    if ($mode =~ /^1/) { &extractSingleKeys(); } else { &extractPairedKeys(); }
}

exit;

###############
# subroutines #
###############
sub usage {
    system("perldoc $0");
    exit;
}

sub loadEnv {
    # loads environment shortcuts contained in the esb-specific shortcuts file
    # essential when running in passive mode in oarsub
    my $f = shift;
    my $perl_command = "perl -MData::Dumper -e 'print Dumper(\\\%ENV)';";
    my $source_line = ". $f 1>&2";
    my %tmp = %{eval('my ' . `$source_line\n$perl_command`)};
    $ENV{$_} = $tmp{$_} for (keys %tmp);
    return 0;
}

sub checkExecutables {
    if (!$cdbyank) { ## make sure the 'cdbyank' executable can be found
	print { *STDERR } BOLD RED 'ERROR:', RESET, qq(\tThe required tool 'cdbyank' cannot be found.\n),
	qq(\tPlease add it to a directory in your path, or provide its exact location with the flag '-c'.\n);
	exit 1;
    }
    if ($bam && !$samtools) {
	print { *STDERR } BOLD RED 'ERROR:', RESET, qq(\tYou indicated the key files are bam alignment files, ),
	qq(but the samtools executables cannot be found.\n\tPlease add them to a directory in your path, or),
	qq(provide their exact location with the flag '-s'.\n);
	exit 1;
    }
    return 0;
}

sub checkInputs {
    ## check for at least one keys file and one index file
    if ( (!$keys1 && !$bam) || !$index1) {
	print { *STDERR } BOLD RED 'ERROR:', RESET,
	qq(\tAt least one Keys file AND one Index file are required.\n);
	exit 1;
    }
    ## make sure that the files provided agree with what the user expects
    my $mode_c = ($keys1 && $keys2) ? ( ($index1 && $index2) ? '22' : '21' ) : 
	( ($index1 && $index2) ? '12' : '11' ); 
    if ($mode) {
	unless ($mode eq $mode_c) {
	    print { *STDERR } BOLD RED 'ERROR:', RESET,
	    qq(\tThe input mode $mode and the mode $mode_c calculated from provided input files do not agree.\n),
	    qq(\tMake sure all intended inputs are included, and that you know your expectations.\n);
	    exit 1;
	}
    } else {
	$mode = $mode_c;
    }
    ## die if asking for bam parsing from two files;
    ## it is not hard to implement for two files, but making it efficient may need to thinking
    if ($mode =~ /^2/ && $bam) {
	print { *STDERR } BOLD RED 'ERROR:', RESET,
	qq(\tExtracting keys from bam files is only available for single key-files.\n);
	exit 1;
    }
    if ($keys1 =~ /\.bam$/ && !$bam) {
	$bam = 1;
	print { *STDERR } BOLD YELLOW 'WARNING:', RESET,
	qq(\tThe 'keys' file carries a '.bam' extension, but the '-b/-bam' flag was not included.\n),
	qq(\tI will assume you meant to include it, and treat your file as a bam file.\n);
    }
    ## check for use of pandaseq-like keys as a single input
    if ($trim_adapt && $mode =~ /^2/ ) {
	print { *STDERR } BOLD RED 'ERROR:', RESET, 
	qq(\tTrimming adaptor sequences from the end of read ids only applies to modes 11 and 12.\n),
	qq(\tPlease revise your entries, or use -h for more information.\n);
	exit 1;
    }
    ## now make sure the input files are actually there
    &checkKeys($keys1);
    &checkKeys($keys2) if ($keys2);
    &checkIndex($index1);
    &checkIndex($index2) if ($index2);

    if ($verbose) {
	my @mode = split(//,$mode);
	my @plu = ( ($mode[0] == 1 ? '' : 's'), ($mode[1] == 1 ? '' : 's') );
	my $act = $trim_adapt ? '' : ' NOT';
	print { *STDERR } BOLD GREEN 'PASS:', RESET, qq(\tAll input files and values seem to be correct.\n),
	qq(\tWill obtain keys from $mode[0] file$plu[0] and extract reads from $mode[1] file$plu[1].\n),
	qq(\tAdaptor indicator sequences will$act be trimmed from the end of read ids.\n),
	qq(\tProceeding with read extractions.\n);
    }
    return 0;
}

sub checkKeys {
    my ($f) = shift;
    unless (-e $f) {
	print { *STDERR } BOLD RED 'ERROR:', RESET, qq(\tCannot find keys in file $f\n);
	exit 1;
    }
    return 0;
}

sub checkIndex {
    my ($f) = shift;
    if ($f =~ /\.cidx$/) {
	if (-e $f) { return 0; }
    } else {
	$f =~ s/\.cidx$//;
	if (-e $f) {
	    system("$ENV{'TOOLS'}/cdbfasta -Q $f");
	    return 0;
	} else {
	    print { *STDERR } BOLD RED 'ERROR:', RESET, qq(\tCannot find either the index file $f.cidx, nor an appropriate reads file\n);
	    exit 1;
	}
    }
}

sub extractPairedKeys {
    my $command = q(awk '!seen[$0]++' ).qq($keys1 $keys2 | $cdbyank);
    &executeExtract($command);
}

sub extractSingleKeys {
    my $command = $bam ? qq($samtools view $keys1 | cut -f 1 | ) : qq(cat $keys1 | );
    $command .= q(sed 's/\:\([A-Z]\+$\)//g' | ) if ($trim_adapt);
    $command .= q(awk '!seen[$0]++' | );
    $command .= $cdbyank;
    &executeExtract($command);
}

sub executeExtract {
    my $cmd = shift;
    my @index = ($index1);
    push(@index, $index2) if ($index2);
    for (my $i=0; $i<@index; $i++) {
	my $out = $keys1;
	my $suf = @index == 1 ? '' : 'R'.($i+1);
	$out =~ s/\.\w+$/\.$suf\.fq/;
	print { *STDERR } BOLD BLACK 'CMD:', RESET, "\t$cmd $index[$i] \> $out\n" if ($verbose);
	system("$cmd $index[$i] \> $out");
    }
}

sub executeExclude {
    my @index = ($index1);
    push(@index, $index2) if ($index2);
    my @reads = ();
    for (my $i=0; $i<@index; $i++) {
	print qq(\t  loading read ids from index file $index[$i]...\n) if ($verbose);
	chomp(my @list = `$cdbyank -l $index[$i]`);
	push(@reads, @list);
    }
    print qq(\t  reducing redundancy in read ids...\n) if ($verbose);
    my %reads = map {$_ => 1} @reads;
    print qq(\t  constructing list of remaining reads...\n) if ($verbose);
    if ($mode =~ /^1/) { %reads = &excludeSingleKeys(\%reads); } else { &excludePairedKeys(\%reads); }

    ## extract and output the remaining reads
    print qq(\t  generating read files with the remaining reads...\n) if ($verbose);
    foreach my $index (@index) {
	my $outfile = $index;
	$outfile =~ s/\.cidx$/\.substracted/;
	my $command = "$cdbyank $index > $outfile";
	open my $cmd, '|-', $command;
	foreach my $k (keys %reads) {
	    print $cmd "$k\n";
	}
    }
}

sub excludeSingleKeys {
    my $hash_ref = shift;

    print BOLD 'PROGRESS:', RESET,
    qq( Excluding reads from a single key file.\n);
    my $command = $bam ? qq($samtools view $keys1 | cut -f 1) : qq(cat $keys1);
    $command .= q( | awk '!seen[$0]++');
    $command .= q( | sed 's/\:\([A-Z]\+$\)//g') if ($trim_adapt);
    open my $cmd, '-|', $command;
    while (<$cmd>) {
	chomp;
	if ($$hash_ref{$_}) { delete($$hash_ref{$_}); } else {
	    print BOLD YELLOW 'WARNING:', RESET,
	    qq( Read $_ from the list of keys was not found in the indexes provided.\n);
	}
    }

    return %{$hash_ref};
}

sub excludePairedKeys {
    my $hash_ref = shift;

    print BOLD 'PROGRESS:', RESET,
    qq( Excluding reads from two key files.\n);
    my $command = q(awk '!seen[$0]++' ).qq($keys1 $keys2);
    open my $cmd, '-|', $command;
    while (<$cmd>) {
	chomp;
	if ($$hash_ref{$_}) { delete($$hash_ref{$_}); } else {
	    print BOLD YELLOW 'WARNING:', RESET,
	    qq( Read $_ from the list of keys was not found in the indexes provided.\n);
	}
    }

    return %{$hash_ref};
}
