#!/usr/bin/perl -w

=head1 NAME

  PandaSeqQC.pl - Parses PANDAseq output, checks for integrity, and gathers statistics

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 COMMENTS

    To run with the LAOTS data, but outside of 02c.preprocess.join_overlaps.sh,
    use the utilities wraper PandaSeqQC.sh; otherwise, required environment
    variables will not be available.

    This script is too clumsy. It works, but it needs to be polished. It is made
    available in its draft form since it is a dependency of 02c.preprocess.join_overlaps.sh
    in the LAO One Year Time Series pipeline.

=head1 AUTHOR - Nicolas Pinel (Institute for Systems Biology - May 30, 2012)

=cut

use strict;
use Getopt::Long;
use POSIX qw(strftime);
use Benchmark qw(timestr timediff) ;

my $max_mismatch_per = 25; ## maximum mismatch percent within overlap to consider the join legitimate

my ($ps_file,$r1_file,$r2_file,$ds,$verbose);
GetOptions('ps=s' => \$ps_file, ## pandaseq output fastq
	   'r1=s' => \$r1_file,   ## index of R1 raw data file, constructed with cdbfasta
	   'r2=s' => \$r2_file,   ## index of R2 raw data file, constructed with cdbfasta
	   'mm=s' => \$max_mismatch_per,
	   'v|verbose' => \$verbose);

&loadENV();

my $date = strftime('%d-%b-%Y %H:%M:%S',localtime); ## date captured for including in the final report

print { *STDERR } qq(Checking for the existence of the input files...\n) if ($verbose);
&checkFile($ps_file);
&checkFile($r1_file);
&checkFile($r2_file);
print { *STDERR } qq(Found all the needed input files...proceeding!\n) if ($verbose);

print { *STDERR } qq(Checking for the presence of index files...) if ($verbose);
my $psix = $ps_file.'.cidx';
my $r1ix = $r1_file.'.cidx';
my $r2ix = $r2_file.'.cidx';
unless ( -e $psix ) { &makeIndex($ps_file); }
unless ( -e $r1ix ) { &makeIndex($r1_file); }
unless ( -e $r2ix ) { &makeIndex($r2_file); }
print { *STDERR } qq(done!\n) if ($verbose);

my %reads = ();
my $t0 = Benchmark->new(); ## see $tf at the end of processing blocs
my $tmp1 = 'tmp1.'.time().'.'.sprintf("%06d",int(rand(999999))); ## the random number prevents accidental file conflicts when
                                                                 ## processing multiple runs within the same directory
my $tmp2 = $tmp1;
$tmp2 =~ s/tmp1/tmp2/;

open my $cmd, '-|', "sed -n \'1~4p\' $ps_file";
open my $t1, '>', $tmp1;  ## contains read headers for extracting reads from PANDAseq fastq
open my $t2, '>', $tmp2;  ## contains read headers for extracting reads from R1 or R2 read files
print { *STDERR } qq(Generating temporary files with read ids...) if ($verbose);
while (<$cmd>) {
  chomp;
  s/^\@//;
  print $t1 $_,"\n";
  s/\:\w+$//;
  print $t2 $_,"\n";
}
close($cmd);
print { *STDERR } qq(done!\n) if ($verbose);

print { *STDERR } qq(Loading reads from the PandaSeq file...) if ($verbose);
my $counter = 0;
open $cmd, '-|', "cat $tmp1 | $ENV{'TOOLS'}/cdbyank $psix";
my $ticker = 0;
my @fields = ();
while (<$cmd>) {

  chomp;
  ++$ticker;
  if ( $ticker == 1 || $ticker == 2 || $ticker == 4 ) {
    push(@fields, $_);
  }

  if ($ticker == 4) {
    my $id = $fields[0];
    $id =~ s/^\@//;
    $id =~ s/\:\w+$//;
    $reads{$id} = {'name'  => $fields[0],
		   'ps_sq' => $fields[1],
		   'ps_qv' => $fields[2]};

    @fields = ();
    $ticker = 0;
  }
}
close($cmd);
undef($cmd);

my $total_reads = scalar keys %reads;
print { *STDERR } qq(loaded $total_reads reads.\n) if ($verbose);

$counter = 0;
$ticker = 0;
@fields = ();
&populateReadPairs($r1ix,1);

$ticker = 0;
$counter = 0;
@fields = ();
&populateReadPairs($r2ix,2);

print { *STDERR } qq(Proceeding with verification of joins...\n\n) if ($verbose);

my $report_file = $ps_file;
$report_file =~ s/\.fastq/\.QC\.report/;
my $misjoin_file = $report_file;
$misjoin_file =~ s/report/misjoined/;
my $fixed_file = $report_file;
$fixed_file =~ s/report/fixed/;

open my $misjoin, '>', $misjoin_file;
open my $fixed, '>', $fixed_file;

my %stats = ('misoverlap' => 0,     ## number of reads flagged as improperly joined
	     'trueoverlap' => 0,    ## number of reads that appear to have been properly joined
	     'miscalls' => 0,       ## mismatches not corrected even though they should have
	     'truecalls' => 0,      ## mismatches that were properly corrected
	     'total_overlap' => 0,  ## total sequence length of overlap across all joined reads
	     'fixed' => 0);         ## reads with miscalls that were fixed (and placed on a separate file)

$counter = 0;
foreach my $k (keys %reads) {
  ++$counter;
  print qq(Verifying join for read $counter of $total_reads\.\n) if ($counter % 100000 == 0);
  
  my $read_length = length($reads{$k}{'r1_sq'});
  my $overlap = 2*$read_length - length($reads{$k}{'ps_sq'});
  $stats{'total_overlap'} += $overlap;
  
  my $ps_s = substr($reads{$k}{'ps_sq'}, ($read_length - $overlap), $overlap);
  my $r1_s = substr($reads{$k}{'r1_sq'}, ($read_length - $overlap), $overlap);
  my $r2_s = substr($reads{$k}{'r2_sq'}, 0, $overlap);
  my $ps_q = substr($reads{$k}{'ps_qv'}, ($read_length - $overlap), $overlap);
  my $r1_q = substr($reads{$k}{'r1_qv'}, ($read_length - $overlap), $overlap);
  my $r2_q = substr($reads{$k}{'r2_qv'}, 0, $overlap);
  
  my %overlaps = ('ps_seq' => [ split(//, $ps_s) ], 
		  'r1_seq' => [ split(//, $r1_s) ],
		  'r2_seq' => [ split(//, $r2_s) ],
		  'ps_qv' => [ split(//, $ps_q) ],
		  'r1_qv' => [ split(//, $r1_q) ],
		  'r2_qv' => [ split(//, $r2_q) ],);
  
  my $mismatches = 0;
  my $miscalls = 0;
  my $truecalls = 0;
  my $uncertaincalls = 0;
  my $fix = 0;
  
  my $alnstr;
  for (my $i=0; $i<$overlap; $i++) {
    
    if ($overlaps{'r1_seq'}[$i] ne $overlaps{'r2_seq'}[$i]) {
      $mismatches++ unless ($overlaps{'r1_seq'}[$i] eq 'N' || $overlaps{'r2_seq'}[$i] eq 'N');
      $alnstr .= ($overlaps{'r1_seq'}[$i] eq 'N' || $overlaps{'r2_seq'}[$i] eq 'N') ? '.' : 'x';
    } else {
      $alnstr .= ' ';
    }
    
  }
  
  my $mismatch_per = sprintf("%.3f", 100*($mismatches/$overlap));
  if ($mismatch_per <= $max_mismatch_per) { 
      $stats{'trueoverlap'}++;
      
      ## fixed miscalled mismatches
      for (my $i=0; $i<$overlap; $i++) {
	  if ($overlaps{'r1_seq'}[$i] eq $overlaps{'r2_seq'}[$i]) {
	      next;
	  } else {
	      my @qvs = (ord($overlaps{'r1_qv'}[$i])-33, ord($overlaps{'r2_qv'}[$i])-33);
	      
	      if ($overlaps{'r1_seq'}[$i] eq $overlaps{'ps_seq'}[$i]) { # PANDAseq chose R1 
		  if ($qvs[0] < $qvs[1]) { ## when it should have picked R2 based on higher QV
		      $miscalls++;
		      $overlaps{'ps_seq'}[$i] = $overlaps{'r2_seq'}[$i];
		      $overlaps{'ps_qv'}[$i] = $overlaps{'r2_qv'}[$i];
		      $fix =1;
		      
		  } elsif ($qvs[0] > $qvs[1]) { ++$truecalls; } else { ++$uncertaincalls; }
		  
	      } elsif ($overlaps{'r2_seq'}[$i] eq $overlaps{'ps_seq'}[$i]) { ## PANDAseq chose R2
		  if ($qvs[0] > $qvs[1]) { ## when it should have picked R1 based on higher QV
		      $miscalls++;
		      $overlaps{'ps_seq'}[$i] = $overlaps{'r1_seq'}[$i];
		      $overlaps{'ps_qv'}[$i] = $overlaps{'r1_qv'}[$i];
		      $fix = 1;
		      
		  } elsif ($qvs[0] < $qvs[1]) { ++$truecalls; } else { ++$uncertaincalls; }
		  
	      }
	      if ($fix) {
		  $stats{'fixed'}++;
		  my $fixed_seq = join('', @{$overlaps{'ps_seq'}});
		  my $fixed_qv = join('', @{$overlaps{'ps_qv'}});
		  $reads{$k}{'ps_sq'} =~ s/$ps_s/$fixed_seq/;
		  
		  ## $reads{$k}{'ps_qv'} =~ s/quotemeta($ps_q)/quotemeta($fixed_qv)/;
		  ## the attempt to adjust the quality value string continues to throw errors
		  ## must find an alternative
		  
		  print $fixed qq($reads{$k}{'name'}\n$reads{$k}{'ps_sq'}\n\+\n$reads{$k}{'ps_qv'}\n);
		  
	      }
	  }
      }
      
  } else {
      $stats{'misoverlap'}++;
    
      my $pad = length($reads{$k}{'ps_sq'}) - length($reads{$k}{'r2_sq'});
      my $str = qq(\>Overlap: $overlap\; Mismatch: $mismatch_per\n);
      $str .= qq($reads{$k}{'name'}\n$reads{$k}{'ps_sq'}\n$reads{$k}{'ps_qv'}\n\n);
      $str .= qq(QVR1\:\t$reads{$k}{'r1_qv'}\nSQR1\:\t$reads{$k}{'r1_sq'}\n);
      $str .= qq(SQPS\:\t$reads{$k}{'ps_sq'}\nSQR2\:\t);
      $str .= ' ' x $pad;
      $str .= qq($reads{$k}{'r2_sq'}\nQVR2\:\t);
      $str .= ' ' x $pad;
      $str .= qq($reads{$k}{'r2_qv'}\n);
      $str .= qq(     \t);
      $str .= ' ' x $pad;
      $str .= qq($alnstr\n\n);
      $str .= '=' x 150;
      $str .= "\n\n";
      print $misjoin $str;
      
  }
  
  $stats{'miscalls'} += $miscalls;
  $stats{'truecalls'} += $truecalls;
  $stats{'uncertain'} += $uncertaincalls;
    
}
print qq(\nFinished all verifications. Now writing report...\n);

my $tf = Benchmark->new();
my $td = timediff($tf,$t0);

my $miscall_frq = sprintf("%.3f", $stats{'miscalls'}/$stats{'total_overlap'});
my $truecall_frq = sprintf("%.3f", $stats{'truecalls'}/$stats{'total_overlap'});

open my $report, '>', $report_file;
print $report "$date\nPandaSeqQC run time: ",timestr($td),"\n";
print $report qq(File:\t$ps_file\nReads:\t$total_reads\n);
print $report qq(Mismatch threshold:\t$max_mismatch_per\n);
print $report qq(Properly joined paired reads\:\t$stats{'trueoverlap'}\n);
print $report qq(Wrongly joined paired reads\:\t$stats{'misoverlap'}\n);
print $report qq(Total overlapping sequence length\:\t$stats{'total_overlap'}\n);
print $report qq(Total repaired mismatches\:\t$stats{'truecalls'}\n);
print $report qq(Total mismatch repair rate\:\t$truecall_frq\n);
print $report qq(Total miscalled mismatches\:\t$stats{'miscalls'}\n);
print $report qq(Total mismatch miscall rate\:\t$miscall_frq\n);
print $report qq(Repaired joined reads\:\t$stats{'fixed'}\n);
print $report qq(Total uncertain mismatches\:\t$stats{'uncertain'}\n\n);

system("rm -f $tmp1 $tmp2");

print qq(All Done!\n\n);

exit;

### subroutines ###

sub usage {

}

sub loadENV {
    my %env = ('TOOLS' => 'directory containing cdbfasta and cdbyank',
	       'LAOTS'  => 'base directory for the One Year Time Series data set');
    
    print { *STDERR } qq(Checking for environment variables...) if ($verbose);
    delete($env{'TOOLS'}) if ($ENV{'TOOLS'} && $ENV{'TOOLS'} ne '');
    delete($env{'LAOTS'}) if ($ENV{'OYTS'} && $ENV{'OYTS'} ne '');
    
    if (scalar keys %env > 0) {
	my $warning = "The following required environment variables were not found:\n";
	foreach my $k (keys %env) {
	    $warning .= "\t$k\t= $env{$k}\n";
	}
	$warning .= "export them to your environment before running this script.\n";
	print { *STDERR } $warning;
	exit 1;
    }
    print { *STDERR } qq(done!\n) if ($verbose);
    return 0;
}

sub checkFile {
    my $f = shift;
    if ( -e $f ) {
	print { *STDERR } qq(  found $f\n) if ($verbose);
	return 0;
    } else {
	print { *STDERR } qq(ERROR: cannot find $f\; aborting the run!\n);
	exit 1;
    }
}

sub makeIndex {
    my $f = shift;
    print { *STDERR } qq(constructing index for $f\...) if ($verbose);
    system("$ENV{'TOOLS'}\/cdbfasta $f -Q");
    print { *STDERR } qq(done!\n) if ($verbose);
}

sub populateReadPairs {
    my $tools = $ENV{'TOOLS'};
    my ($file, $pair) = @_;
    
    print { *STDERR } qq(Loading reads from $file\...) if ($verbose); 
    
    open my $fh, '-|', "cat $tmp2 | $tools/cdbyank $file";
    while (<$fh>) {
	chomp;
	++$ticker;
	
	if ( $ticker == 1 || $ticker == 2 || $ticker == 4 ) {
	    
	    if ($ticker == 1) {
		++$counter;
		push(@fields, $_);
		print qq(Finding in R$pair: read $counter of $total_reads\.\.\.\n) if ($counter % 100000 == 0);
		
	    } elsif ($ticker == 2) { ## reverse/complement
		if ($pair == 2) {
		    my @str = split(//);
		    @str = reverse(@str);
		    $_ = join('',@str);
		    tr/ACGT/TGCA/;
		}
		push(@fields, $_);
		
	    } elsif ($ticker == 4) {
		if ($pair == 2) {
		    my @str = split(//);
		    @str = reverse(@str);
		    $_ = join('',@str);
		}
		push(@fields, $_);
		
		my $id = $fields[0];
		$id =~ s/^\@//;
		$id =~ s/ \d{1}\:.*//;
		$reads{$id}{"r$pair\_sq"} = $fields[1];
		$reads{$id}{"r$pair\_qv"} = $fields[2];
		
		@fields = ();
		$ticker = 0;
	    }
	    
	}
    }
    print qq(Finished loading the R$pair pair for $counter reads.\n\n);
    
    return 0;
}
