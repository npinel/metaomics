#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 25 October 2013                                                                  # 
# Script: 05a.read_mapping.bwa.sh                                                        #
# Purpose:                                                                               #
#    Bash shell wrapper script to configure and run bwa alignments                       #
# NOTES:                                                                                 #
#    Full description of the command line options is available from:                     #
#    http://bio-bwa.sourceforge.net/bwa.shtml                                            #
##########################################################################################

#
# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldblk='\e[1;30m' # bold black
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

# load ESB shortcuts
if [[ -z $ESB ]]; then
    source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
    export PATH=$PATH:$TOOLS
fi

# hard code path to desired versions of bwa/samtools
# should use modules, but have not figured out how to make them work
BWA=$TOOLS"/bwa-0.7.8/bin/bwa"
SAMTOOLS=$TOOLS"/samtools-0.1.19/bin/samtools"
#BWA=`which bwa`
#SAMTOOLS=`which samtools`

if [[ -z $BWA ]]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}An executable for bwa cannot be found."
    exit 1
fi

if [[ -z $SAMTOOLS ]]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}An executable for samtools cannot be found."
    exit 1
fi

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
  -m  mapping quality threshold for filtering output bam file
  -u  also create bam file with unmapped reads, for subsequent extraction
  -x  exclude alignment bam file (use in conjunction with -u if only unmapped reads are desired)
  -s  include single-end reads
  -p  include pandaseq reads in the single-end data set
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

build_index(){
echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Constructing an index for the reference sequences."
    $BWA index $REFSEQ
}

#create_log(){

#cmd='/mnt/nfs/projects/ecosystem_biology/local_tools/bwa --version'
#LOG=$WORKDIR/$0.`whoami`.`date +%Y%m%d_%T`.log
#`cp $0 $LOG`
#`$cmd >> $LOG`
#`chmod 640 $LOG`
#    echo 'creating log'
#}

filter_sam(){
    out=$1
    echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Filtering output SAM file to exclude unmapped reads."
    echo `$SAMTOOLS view -bhS -F 0x0004 $out.sam | samtools sort - $out`
    wait ${!}
    echo `rm $out.sam`
}

####################
# define variables #
####################
FILE=
DATA_TYPE=MG
REFSEQ=
REFDESC=
QUAL=
ODIRSUFF=
INCLUDE_PS=0
INCLUDE_SE=0
OUTPUT_ALIGNED=1
OUTPUT_UNMAPPED=0
WORKDIR=
TESTING=0
MAPQ=0

while getopts "f:r:i:q:d:w:m:pushxt" OPTION
do
    case $OPTION in
	h) usage ;;
	f) FILE=$OPTARG ;;
	r) REFSEQ=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	i) REFDESC=$OPTARG ;;
	m) MAPQ=$OPTARG ;;
	u) OUTPUT_UNMAPPED=1 ;;
	x) OUTPUT_ALIGNED=0 ;;
	p) INCLUDE_PS=1 ;;
	s) INCLUDE_SE=1 ;;
	d) DATA_TYPE=$OPTARG ;;
	w) WORKDIR=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

if [[ -z "$FILE" ]]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}A sample list file, or a sample identifier are required."
    usage
fi

THREADS=`cat $OAR_RESOURCE_PROPERTIES_FILE | wc -l` # the file contains 1 line per core 
if [[ -z $THREADS ]]; then THREADS=1; fi

# make sure reference sequence has been indexed
if [ ! -f "$REFSEQ.amb" ]; then build_index $REFSEQ; fi

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
    DATDIR=$LAOTS/$S/$DATA_TYPE
    WORKDIR=/work/users/`whoami`/$S.bwa.$WORKDIR
    if [[ ! -d $WORKDIR ]]; then
	mkdir $WORKDIR
    fi

    R1=$DATDIR/${S}_$DATA_TYPE.$QUAL.R1.fq
    R2=$DATDIR/${S}_$DATA_TYPE.$QUAL.R2.fq
    SE=$DATDIR/${S}_$DATA_TYPE.$QUAL.SE.fq

    if [ -f $R1 ] && [ -f $R2 ]; then
	SAMHEADER="@RG\tID:$S.$DATA_TYPE.$QUAL\tLB:$S.$DATA_TYPE.$QUAL\tSM:$S.$DATA_TYPE.$QUAL\tPL:Illumina\tCN:TGen"

	# 1.) perform the mapping
	echo -e "${bldblk}[`date +%H:%M:%S` ${bldgrn}START${bldblk}] ${txtrst}Processing ${bldblk}$S${txtrst}."
	echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Indexing and aligning paired reads..."
	$BWA mem -v 1 -t $THREADS -M -R $SAMHEADER $REFSEQ  $R1 $R2 > $WORKDIR/${S}_$DATA_TYPE.bwa.aln.pairs.sam

	# now work on single end reads if instructed to do so
	if [[ $INCLUDE_SE -eq 1 ]]; then
	    echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Indexing and aligning single-end reads..."
	    $BWA mem -v 1 -t $THREADS -M -R $SAMHEADER $REFSEQ  $SE > $WORKDIR/${S}_$DATA_TYPE.bwa.aln.se.sam
	fi

	# 2.) output alignment files
	if [[ $OUTPUT_ALIGNED -eq 1 ]]; then
	    echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Sorting SAM file for paired reads..."
	    $SAMTOOLS view -ubhS -F 0x0004 -q $MAPQ $WORKDIR/${S}_$DATA_TYPE.bwa.aln.pairs.sam | samtools sort - $WORKDIR/${S}_$DATA_TYPE.bwa.aln.pairs
	    if [[ $INCLUDE_SE -eq 1 ]]; then
		echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Sorting SAM file for single-end reads..."
		$SAMTOOLS view -ubhS -F 0x0004 -q $MAPQ $WORKDIR/${S}_$DATA_TYPE.bwa.aln.se.sam | samtools sort - $WORKDIR/${S}_$DATA_TYPE.bwa.aln.se
		echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Merging BAM files..."
		OUTBAM=$DATDIR/${S}_${DATA_TYPE}.bwa.${REFDESC}_x_$QUAL.mapq$MAPQ.bam
		$SAMTOOLS merge $OUTBAM $WORKDIR/${S}_$DATA_TYPE.bwa.aln.*.bam
	    else
		OUTBAM=$DATDIR/${S}_${DATA_TYPE}.bwa.${REFDESC}_x_$QUAL.mapq$MAPQ.pairs_only.bam
		mv $WORKDIR/${S}_$DATA_TYPE.bwa.aln.pairs.bam $OUTBAM
	    fi
	fi

	if [[ $OUTPUT_UNMAPPED -eq 1 ]]; then
	    echo -e "${bldblk}[`date +%H:%M:%S`  ${bldblu}PROG${bldblk}] ${txtrst}Identifying unmapped reads..."
	    $SAMTOOLS view -bhS -f 0x0004 $WORKDIR/${S}_$DATA_TYPE.bwa.aln.pairs.sam > $DATDIR/${S}_${DATA_TYPE}.bwa.${REFDESC}_x_$QUAL.unmapped.bam
	fi

	rm $WORKDIR/${S}_$DATA_TYPE.bwa.aln.*.*am
	rmdir $WORKDIR
	echo -e "${bldblk}[`date +%H:%M:%S`   ${bldgrn}END${bldblk}] ${txtrst}Completed ${bldblk}$S${txtrst}."
    fi
done

exit 0
