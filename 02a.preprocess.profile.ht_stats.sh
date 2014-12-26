#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 8 April 2013                                                                     #
# Script: 02a.preprocess.profile.ht_stats.sh                                             #
# Purpose:                                                                               #
#    Batch profiling of fastq-formatted sequence reads, in linear execution, performed   #
#    with ht_stats
#                                                                                        #
##########################################################################################

# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldblk='\e[1;30m' # bold black
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

####################
# DEFINE FUNCTIONS #
####################
usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -d  data type; values = MG|MT (required; default MG)
  -f  file with sample list, or single dataset name (required)
  -l  maximum read length (default 105 nt)
  -q  qualifier for the read set (e.g., if a read file is D01_MG.clp7.uniq.R1.fq, set -q clp7.uniq)
  -p  assume paired input files, flag
  -t  test, flag
  -h  display this message

EOF
    exit 1
}

rename_output(){
    outdir=$1
    sample=$2
    type=$3

    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Renaming output files found in $outdir."
    OUTFILES=$outdir/*
    for of in $OUTFILES
    do
	nf=$outdir/${sample}_$type.`basename $of`
	`mv $of $nf`
    done
}

################
# LOAD MODULES #
################

## no modules have yet been greated for HTQC, so run from the $TOOLS directory

####################
# DEFINE VARIABLES #
###################
# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
export PATH=$PATH:$TOOLS

FILE=
TYPE=MG
PAIRED=
QUAL=
RLENGTH=105
TESTING=

while getopts "d:f:l:q:pth" OPTION
do
    case $OPTION in
	d) TYPE=$OPTARG ;;
	f) FILE=$OPTARG ;;
	l) RLENGTH=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	p) PAIRED=1 ;;
	t) TESTING=1 ;;
	?) usage ;;
	h) usage ;;
    esac
done

if [[ -z $FILE ]] || [[ -z $TYPE ]]
then
    usage
fi

if [[ -n $QUAL ]]; then
    QUAL='.'$QUAL
fi

THREADS=`cat $OAR_RESOURCE_PROPERTIES_FILE | wc -l` # the file contains 1 line per core
if [ -z $THREADS ]; then THREADS=1; fi

#######################
# identify input type #
#######################

# declare the array
declare -a SAMPLES=()
if [ -f $FILE ] ; then
    while read -r line; do
	[[ "$line" =~ ^#.*$ ]] && continue
	set -- $line
	declare -a fields=($*)
	SAMPLES=("${SAMPLES[@]}" ${fields[0]})
    done < $FILE
else
    SAMPLES=($FILE)
fi

#####################
# START THE PROCESS #
#####################
for S in ${SAMPLES[@]} ; do # loop through array of sample id's
    echo -e "${bldblk}[`date +%H:%M:%S`    ${bldgrn}START${bldblk}] ${txtrst}Starting processing reads for ${bldblk}$S${txtrst}."
    R1=$LAOTS/$S/$TYPE/${S}_${TYPE}$QUAL.R1.fq
    R2=$LAOTS/$S/$TYPE/${S}_${TYPE}$QUAL.R2.fq

    if [[ -n $R1 ]] && [[ -n $R2 ]] ; then
	echo -e "> Found both R1 and R2 files for the indicated reads."
	if [[ -n $PAIRED ]] ; then
	    OUTDIR=$LAOTS/$S/$TYPE/qc.ht_stat$QUAL.pairs
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Profiling read pairs."
	    echo -e "> $TOOLS/ht_stat --encode=sanger -q -p -l $RLENGTH -t $THREADS -m F -o $OUTDIR $R1 $R2"
	    `$TOOLS/ht_stat --encode=sanger -q -p -l $RLENGTH -t $THREADS -m F -o $OUTDIR $R1 $R2`
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Finished collecting statistics. Now constructing plots."
	    `$TOOLS/ht_stat_draw.pl --dir $OUTDIR`
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Finished plotting."
	    rename_output $OUTDIR $S $TYPE	    
	else 
	    echo -e "> However, the -p flag was not set, so will profile read files individually."
	    # process R1
	    OUTDIR=$LAOTS/$S/$TYPE/qc.ht_stat$QUAL.R1
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Profiling reads in R1"
	    echo -e "> $TOOLS/ht_stat --encode=sanger -q -l $RLENGTH -t $THREADS -m F -o $OUTDIR $R1"
	    `$TOOLS/ht_stat --encode=sanger -q -l $RLENGTH -t $THREADS -m F -o $OUTDIR $R1`
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Finished collecting statistics. Now constructing plots."
	    `$TOOLS/ht_stat_draw.pl --dir $OUTDIR`
	    rename_output $OUTDIR $S $TYPE
	    # process R2
	    OUTDIR=$LAOTS/$S/$TYPE/qc.ht_stat$QUAL.R2
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Profiling reads in R2"
	    echo -e "> $TOOLS/ht_stat --encode=sanger -q -l $RLENGTH -t 12 -m F -o $OUTDIR $R2"
	    `$TOOLS/ht_stat --encode=sanger -q -l $RLENGTH -t 12 -m F -o $OUTDIR $R2`
	    echo -e "${bldblk}[`date +%H:%M:%S` PROGRESS] ${txtrst}Finished collecting statistics. Now constructing plots."
	    `$TOOLS/ht_stat_draw.pl --dir $OUTDIR`
	    rename_output $OUTDIR $S $TYPE
	fi
    else
	echo -e "${bldblk}[`date +%H:%M:%S`  ${bldred}WARNING${bldblk}] ${txtrst}Could not find the read files for ${bldblk}$S${txtrst}. Skipping."
    fi
    echo -e "${bldblk}[`date +%H:%M:%S`      ${bldgrn}END${bldblk}] ${txtrst}Finished processing reads for ${bldblk}$S${txtrst}."
done

##########
# notify #
##########
#DATE=`date`
#`mail -r "npinel@gaia.uni.lu" -s "test" "pinel.n@gmail.com" << END_MAIL
#All the files have been processed by ht_stat
#Terminating job, $DATE
#END_MAIL`

exit
