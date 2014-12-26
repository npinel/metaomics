#!/bin/bash

##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 22 April 2013                                                                    #
# Script: PandaSeqQC.sh                                                                  #
# Purpose:                                                                               #
#    To run the Panda Seq Quality Control script (PandaSeqQC.pl) with OYTS data          #
#    but outside of the 02c.preprocess.join_overlaps.sh                                  #
#    Needed in order to feed the correct environment variables to the perl script        #
# NOTES:                                                                                 #
##########################################################################################

usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -d  data set from the One Year Time Series (required)
  -t  data type; values = Metagenomic|MG|Metatranscriptomic|MT (required)
  -m  maximum mismatch in overlapping region (default = 25)
  -v  display run time information (verbose)
  -h  display this message
    
EOF
    exit 1
}

DATA_SET= ;DATA_TYPE= ; MAX_MISMATCH= ;VERBOSE=

while getopts "d:t:m:vh" OPTION
do
    case $OPTION in
	h) usage ;;
	d) DATA_SET=$OPTARG ;;
	t) DATA_TYPE=$OPTARG ;;
	m) MAX_MISMATCH=$OPTARG ;;
	v) VERBOSE=1 ;;
	?) usage ;;
    esac
done

if [[ -z $DATA_SET ]] || [[ -z $DATA_TYPE ]]
then
    usage
fi

source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

PS=$LAOTS/$DATA_SET/$DATA_TYPE/${DATA_SET}_${DATA_TYPE}_PS.raw.fastq
R1=$LAOTS/$DATA_SET/$DATA_TYPE/${DATA_SET}_${DATA_TYPE}_R1.clp7.fastq
R2=$LAOTS/$DATA_SET/$DATA_TYPE/${DATA_SET}_${DATA_TYPE}_R2.clp7.fastq

COMMAND="-ps $PS -r1 $R1 -r2 $R2"
if [[ -n $MAX_MISMATCH ]]; then COMMAND=$COMMAND" -mm $MAX_MISMATCH"; fi
if [[ -n $VERBOSE ]]; then COMMAND=$COMMAND" -v"; fi

$TOOLS/PandaSeqQC.pl $COMMAND
