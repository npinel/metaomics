#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 3 November 2013                                                                  # 
# Script: 05a.read_mapping.bowtie2.sh                                                    #
# Purpose:                                                                               #
#    Bash shell wrapper script to configure and run bowtie2 alignments                   #
# NOTES:                                                                                 #
#    Full description of the command line options is available from:                     #
#    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml                              #
##########################################################################################

# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

# define functions
usage(){
    cat << EOF

Usage: $0 options

OPTIONS:
  -f  file with tab-separated sample list (ID\tYYYY-MM-DD; required)
  -r  reference sequence file (full path; required)
  -i  reference description (for appending to output sam/bam file name)
  -d  data type (MG vs MT; default=MG)
  -q  qualifier string for read files (e.g., if read file is D11_MG.clp7.no_ovlp.qv20.lng40.R1.fq,
      the qualifier corresponds to 'clp7.no_ovlp.qv20.lng40')
  -c  comment to be included in the log file
  -p  include pandaseq reads in the single-end data set
  -s  direct bowtie2 to output mapped/unmapped reads in new fastq files
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

build_index(){
    bwtbld=`which bowtie2-build`
    if [ -z "$bwtbld" ]; then
	echo -e "${bldred}ERROR${txtrst}:\tA bowtie2-formatted index for the reference is missing."
	echo -e "\tCannot construct it because an executable for bowtie2-build cannot be found."
	exit 1
    else
	refseq=$1
	echo `$bwtbld $refseq $refseq 1>&2 /dev/null`
    fi
}

create_log(){

#cmd='/mnt/nfs/projects/ecosystem_biology/local_tools/bowtie2 --version'
#LOG=$WORKDIR/$0.`whoami`.`date +%Y%m%d_%T`.log
#`cp $0 $LOG`
#`$cmd >> $LOG`
#`chmod 640 $LOG`
    echo 'creating log'
}

filter_sam(){
    out=$1
    echo "filtering output SAM file to exclude unmapped reads"
    echo `samtools view -bhS -F 0x0004 $out.sam | samtools sort - $out`
    wait ${!}
    echo `rm $out.sam`
}

####################
# define variables #
####################
FILE=
DATA_TYPE=MG
THREADS=12
REF=
REFDESC='ref'
INCLUDE_PS=0
TESTING=0
SAVE_READS=0
while getopts "f:r:i:q:d:c:psht" OPTION
do
    case $OPTION in
	h) usage ;;
	f) FILE=$OPTARG ;;
	r) REF=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	i) REFDESC=$OPTARG ;;
	p) INCLUDE_PS=1 ;;
	c) COMMENT=$OPTARG;;
	d) DATA_TYPE=$OPTARG ;;
	s) SAVE_READS=1 ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
export PATH=$PATH:$TOOLS

module load Bowtie2/2.0.2-ictce-5.3.0

bwt=`which bowtie2`
if [ -z "$bwt" ]; then
    echo -e "${bldred}ERROR${txtrst}:\tAn executable for bowtie2 cannot be found."
    exit 1
fi

if [[ -z $FILE ]]; then
    echo -e "${bldred}ERROR${txtrst}:\tA sample list file, or a sample identifier are required."
    usage
fi

# make sure reference sequence has been indexed
if [ ! -f "$REF.1.bt2" ]; then build_index $REF; fi

#################################################
# identify input type and populate sample array #
#################################################
declare -a SAMPLES=()
if [ -f $FILE ]
then
    while read -r line; do
	[[ "$line" =~ ^#.*$ ]] && continue
	set -- $line
	declare -a fields=($*)
	SAMPLES=("${SAMPLES[@]}" ${fields[0]})
    done < $FILE
else
    SAMPLES=("${SAMPLES[@]}" $FILE)
fi

#######################
# execute the mapping #
#######################
for S in "${SAMPLES[@]}"
do

    WORKDIR=$LAOTS/$S/$DATA_TYPE
    DATDIR=$WORKDIR
    SAM=$WORKDIR/${S}_${DATA_TYPE}.mapped_$REFDESC # un-sorted sam file output
    
    # paired end reads (as parallel fastq files)
    R1=$DATDIR/${S}_$DATA_TYPE.$QUAL.R1.fq
    R2=$DATDIR/${S}_$DATA_TYPE.$QUAL.R2.fq

    if [ -e $R1 ] && [ -e $R2 ] ; then
        # comma-separated list of single-end reads (e.g., PandaSeq joins, singletons)
	SE=$DATDIR/${S}_$DATA_TYPE.$QUAL.SE.fq
	if [[ $INCLUDE_PS -eq 1 ]]; then
	    SE=$SE,$DATDIR/${S}_${DATA_TYPE}_PS.fixed.fastq ## FIX: this is a very ugly hack
	fi

	cmd=$bwt
	cmd=$cmd' --very-sensitive-local' # preset mode equivalent to: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
	cmd=$cmd" -q -t -p $THREADS"      # -q to indicate reads are in fastq format
                                          # -t to output time information
                                          # -p 12 for 12 cores (1 node in gaia)
	cmd=$cmd" -x $REF"                # base of the index files, without the '.1.bt2,.rev.1.bt2' extensions
	cmd=$cmd" -1 $R1 -2 $R2"          # paired-end reads, in separate (not interleaved) files
	cmd=$cmd" -U $SE"                 # files containing single-end reads
	cmd=$cmd" -S $SAM.sam"            # output (unsorted) sam file
	cmd=$cmd" --rg-id $S.$DATA_TYPE.$QUAL"
	if [[ $SAVE_READS -eq 1 ]]; then
        # if indicated, directing bowtie2 to save
        # mapped/unmapped reads
	    cmd=$cmd' --un $WORKDIR/${S}_${DATA_TYPE}_x_$REFDESC.btw2.SE.Unmapped.fastq'
	    cmd=$cmd" --al $WORKDIR/${S}_${DATA_TYPE}_x_$REFDESC.btw2.SE.Mapped.fastq"
	    cmd=$cmd" --un-conc $WORKDIR/${S}_${DATA_TYPE}_x_$REFDESC.btw2.PE.Disconcordant_R%.fastq"
	    cmd=$cmd" --al-conc $WORKDIR/${S}_${DATA_TYPE}_x_$REFDESC.btw2.PE.Concordant_R%.fastq"
	fi

	create_log $cmd # create a log of the run values    
	echo $cmd
	$cmd

        # filtering the sam file to retain only mapped reads, and store as bam format
	filter_sam $SAM
    fi
done

exit 0
