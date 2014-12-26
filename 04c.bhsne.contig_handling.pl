#!/usr/bin/perl -w

=head1 NAME

  04c.bhsne.contig_handling.pl - prepare contigs for BH-SNE, or extract from indeces

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHOR - Nicolas Pinel (29 April 2013)

=cut

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Bio::SeqIO;

my ($seq_file,$dir_file,$work_dir,$verbose,$help);
GetOptions('i|input=s'    => \$seq_file,
	   'd|seq_dir=s'  => \$dir_file,
	   'w|work_dir=s' => \$work_dir,
	   'v|verbose'    => \$verbose,
	   'h|help'       => \$help);

chomp(my $cdbfasta = `which cdbfasta`);
chomp(my $cdbyank = `which cdbyank`);
chomp(my $pwd = `pwd`);

my $sequences = &buildSeqHash();

if ($seq_file && $dir_file) {
    &extractSequences($sequences);
} else {
    &constructSeqDir($sequences);
}

exit;

###############
# subroutines #
###############
sub buildSeqHash {
    my %hash = ();
    my $cnt = 0;
    my $seqin = Bio::SeqIO->new(-file => $seq_file,
				-format => 'fasta');
    while (my $seqobj = $seqin->next_seq()) {
	$hash{++$cnt} = $seqobj->display_id();
    }
    return \%hash;
}

sub extractSequences {
    my $seq_hash = shift;
    $work_dir = $work_dir ? $work_dir : $pwd;
    my @subsets = <$work_dir/Subset*\_indices.txt>;
    unless (-e "$seq_file\.cidx") { system("$cdbfasta $seq_file"); } ## maybe the index file should be on work_dir?

    foreach my $set (@subsets) {
	my $ofile = $set;
	$ofile =~ s/\_indices.txt$/\_sequences.fa/;
	my $command = "$cdbyank $seq_file\.cidx > $ofile";
	open my $cmd, '|-', $command;
	open my $in, '<', $set;
	while (<$in>) {
	    chomp;
	    print $cmd "$$seq_hash{$_}\n";
	}
    }
}

sub constructSeqDir {
    my $seq_hash = shift;
    chomp(my $ofile = `basename $seq_file`);
    $ofile =~ s/\.\w+$/\.seqdir/;
    $ofile = $pwd.'/'.$ofile;
    open my $o, '>', $ofile;
    foreach my $k (keys %{$seq_hash}) {
	print $o qq($$seq_hash{$k}\t$k\n);
    }
    return 0;
}
