#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 8 April 2013                                                                     #
# Script: 02c.preprocess.join_overlaps.sh                                                #
# Purpose:                                                                               #
#    Join overlapping members of a read pair using PANDASeq (PS), and apply quality      #
#    corrections to the output. It 're-splits' the joined reads in order to use them     #
#    with assemblers that require paired ends reads (e.g., idba_ud)                      #
# NOTES:                                                                                 #
#    As of Apr 15, it is still missing the last three steps:
#      1. fix pandaseq output files
#      2. extract joined reads from the R1/R2 files
#      3. split joined reads to create pseudo-paired end files
##########################################################################################

# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -t  data type; values = Metagenomic|MG|Metatranscriptomic|MT (required)
  -f  file with sample list, or single dataset name (required)
  -h  display this message

EOF
    exit 1
}

run_pandaseq(){
# function: run_pandaseq
# purpose: to perform pair joining using PANDASeq
    testing=$1
    if [[ $testing -eq 1 ]]
    then
	echo "testing: $TOOLS/pandaseq -f $2 -r $3 -o 25 -t 0.9 -F 1> $4 2> $5"
    else
	echo "executing: $TOOLS/pandaseq -f $2 -r $3 -o 25 -t 0.9 -F 1> $4 2> $5"
	`$TOOLS/pandaseq -f $2 -r $3 -o 25 -t 0.9 -F 1> $4 2> $5`
    fi
}

index_fastq(){
# function: index_fastq
# purpose: to index fastq files so that the QC script
# may extract reads from indexed files rapidly;
# it will only run if it fails to detect the index files,
# the default for which appends '.cidx' to the fastq file.
    testing=$1
    file=$2
    if [ -e $file.cidx ]
	then
	echo "Index exists for $file"
	echo "Moving on..."
    else
	if [[ $testing -eq 1 ]]
	then
	    echo "testing: $TOOLS/cdbfastq -Q $file"
	else
	    echo "executing: $TOOLS/cdbfasta -Q $file"
	    `$TOOLS/cdbfasta -Q $file`
	fi
    fi
}

qc_pandaseq(){
# function: qc_pandaseq
# purpose: runs custom script PandaSeqQC.pl to conduct a quality control
# of joins performed by PANDASeq
    testing=$1
    index_fastq $testing $2
    index_fastq $testing $3
    index_fastq $testing $4
    
    if [[ $testing -eq 1 ]]
    then
	echo "testing: $TOOLS/PandaSeqQC.pl -ps $4 -psix $4.cidx -r1ix $2.cidx -r2ix $3.cidx -mm 25"
    else
	echo "executing: $TOOLS/PandaSeqQC.pl -ps $4 -psix $4.cidx -r1ix $2.cidx -r2ix $3.cidx -mm 25"
	`$TOOLS/PandaSeqQC.pl -ps $4 -psix $4.cidx -r1ix $2.cidx -r2ix $3.cidx -mm 25`
    fi
}

####################
# DEFINE VARIABLES #
####################
PWD=`pwd`
FILE=
DATA_TYPE=MG
TESTING=0
while getopts "f:d:ht" OPTION
do
    case $OPTION in
	h) usage ;;
	d) DATA_TYPE=$OPTARG ;;
	f) FILE=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

if [[ -z $FILE ]]
then
    usage
fi

# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

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

####################################
# loop through samples and process #
####################################
for S in "${SAMPLES[@]}"
do
    WD=$PWD/$S/$DATA_TYPE
    OFILE=$WD/${S}_${DATA_TYPE}_RawPS.fastq
    LOG=$WD/${S}_${DATA_TYPE}_RawPS.log
    R1=$WD/${S}_${DATA_TYPE}_R1.clp7.fastq
    R2=$WD/${S}_${DATA_TYPE}_R2.clp7.fastq
    
    run_pandaseq $TESTING $R1 $R2 $OFILE $LOG
    qc_pandaseq $TESTING $R1 $R2 $OFILE
done

exit 0
