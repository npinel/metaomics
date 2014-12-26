#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 27 November 2013                                                                 # 
# Script: util.count_reads.sh                                                            #
# Purpose:                                                                               #
#    Bash shell scripts to count reads in fq files in batch                              #
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

# define functions
usage(){
    cat << EOF

Usage: $0 options

OPTIONS:
  -f  file with tab-separated sample list (ID\tYYYY-MM-DD; required)
  -d  data type (MG vs MT; default=MG)
  -q  qualifier string for read files (e.g., if read file is D11_MG.clp7.no_ovlp.qv20.lng40.R1.fq,
      the qualifier corresponds to 'clp7.no_ovlp.qv20.lng40')
  -s  include single-end reads
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

####################
# define variables #
####################
FILE=
DATA_TYPE=MG
QUAL=
INCLUDE_SE=0
TESTING=0
MAPQ=5
while getopts "f:q:d:sht" OPTION
do
    case $OPTION in
	f) FILE=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	s) INCLUDE_SE=1 ;;
	d) DATA_TYPE=$OPTARG ;;
	t) TESTING=1 ;;
	h) usage ;;
	?) usage ;;
    esac
done

# load ESB shortcuts
if [[ -z $ESB ]]; then
    source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
fi

if [ -z $FILE ]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}A sample list file, or a sample identifier are required."
    usage
fi

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

########################
# execute the counting #
########################
for S in "${SAMPLES[@]}"
do
    DATDIR=$LAOTS/$S/$DATA_TYPE
    R1=$DATDIR/${S}_$DATA_TYPE.$QUAL.R1.fq
    R2=$DATDIR/${S}_$DATA_TYPE.$QUAL.R2.fq
    SE=$DATDIR/${S}_$DATA_TYPE.$QUAL.SE.fq

    if [ -f $R1 ] && [ -f $R2 ]; then
	let "WC_R1=`wc -l $R1 | cut -f 1 -d ' '`/2"
#	let "WC_R2=`wc -l $R2 | cut -f 1 -d ' '`/4"
	# now work on single end reads if instructed to do so
	if [[ $INCLUDE_SE -eq 1 ]]; then
	    let "WC_SE=`wc -l $SE | cut -f 1 -d ' '`/4"
	    echo -e "$S\t$WC_R1\t$WC_SE"
	else
	    echo -e "$S\t$WC_R1"
	fi
    fi
done

exit 0
