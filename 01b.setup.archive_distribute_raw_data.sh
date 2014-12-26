#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 9 April 2013 (last edited on 18 November 2013)                                   #
# Script: 01b.setup.archive_dristribute_raw_data.sh                                      #
# Purpose:                                                                               #
#    Decompress, distribute, and re-compress/archive raw data from TGen-delivered system #
#    to the Gaia $LAOTS file system.                                                     #
# NOTES:                                                                                 #
#    1. This script assumes a standard notation for the data files as delivered by TGen. #
#       Such format must include in the name a 'D##' string that identifies the sample   #
#       number (e.g., D2, D35, etc); and a R1 or R2 string that identifies the           #
#       orientation of the reads.                                                        #
#    2. The script will also create a symbolic link within each of the sample work       #
#       directories of the form D##_MX.R#.fq, which is the naming convention for down-   #
#       stream scripts in the pipeline.                                                  #
#    3. A log file of name log.01b.setup.archive_distribute_raw_data.sh.<date>_<time>    #
#       will be created in the directory whence this script is invoked.                  #
##########################################################################################

# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldbck='\e[1;30m' # bold black
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -d  directory where *fastq.gz files to be transferred are located (required)
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

get_sample() {
    f=$1
    sample=
    if [[ $f =~ D[0-9]{1,2} ]]; then
	sample=${BASH_REMATCH[0]}
	if [[ $sample =~ [0-9]+ ]]; then
	    number=$(printf "%02d" ${BASH_REMATCH[0]})
	    sample=D$number
	fi
    else
	echo "A standard sample code/id was not identified in $f." 1>&2
	exit 1
    fi
    echo $sample
}

get_orientation() {
    f=$1
    r=
    if [[ $f =~ R[0-9]{1} ]]; then
	r=${BASH_REMATCH[0]}
    else
	echo "A standard orientation notation (R1 or R2) was not identified in $f." 1>&2
	exit 1
    fi
    echo $r
}

get_type_long() {
    f=$1
    data=
    if [[ $f =~ 'Metagenom' ]]; then
	data='Metagenomic'
    elif [[ $f =~ 'Metatranscriptom' ]]; then
	data='Metatranscriptomic'
    else
	echo "The type of data (Metagenomic or Metatranscriptomic) was not understood from the file name in $f." 1>&2
	exit 1
    fi
    echo $data
}

get_type_short() {
    d=$1
    t=
    if [ "$d" = 'Metagenomic' ]; then t='MG'; else t='MT'; fi
    echo $t
}

####################
# DEFINE VARIABLES #
####################
DIR= # directory containing the data to be distributed and archived
TESTING=
while getopts "d:th" OPTION
do
    case $OPTION in
	d) DIR=$OPTARG ;;
	t) TESTING=1 ;;
	h) usage ;;
	?) usage ;;
    esac
done

if [[ -z $DIR ]]; then usage ; fi

source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

log_file=`pwd`/log.`basename $0`.`date +%Y%m%d_%H%M%S`
dir_scratch=$SCRATCH/lao_archive_distribute_raw_data.`date +%s`
mkdir $dir_scratch

##############################
# gather a list of the files #
##############################
declare -a FILES=$DIR/*.fastq.gz
for file in ${FILES[@]}; do
    SAMPLE=$(get_sample $file)
    R=$(get_orientation $file)
    DATA=$(get_type_long $file)
    TYPE=$(get_type_short $DATA)

    echo -e "cp $file $dir_scratch" >> $log_file
    if [ -z $TESTING ]; then cp $file $dir_scratch ; fi
    file_scratch=$dir_scratch/$(basename $file)
    echo -e "gzip -d $file_scratch" >> $log_file
    if [ -z $TESTING ]; then gzip -d $file_scratch ; fi
    file_scratch=${file_scratch/\.gz/}
    dir_work=$LAOTS/$SAMPLE/$TYPE
    echo -e "cp $file_scratch $dir_work/" >> $log_file
    if [ -z $TESTING ]; then cp $file_scratch $dir_work/ ; fi
    file_work=$(basename $file_scratch)
    symlink=${SAMPLE}_${TYPE}.$R.fq
    echo -e "cd $dir_work" >> $log_file
    cd $dir_work
    echo -e "ln -s $file_work $symlink" >> $log_file
    if [ -z $TESTING ]; then ln -s $file_work $symlink ; fi
    echo >> $log_file
done

################################
# decompress files and archive #
################################
declare -a FILES=$dir_scratch/*.fastq
for file in ${FILES[@]}; do
    SAMPLE=$(get_sample $file)
    DATA=$(get_type_long $file)
    dir_arch=$ARCHIVE/LAO/$DATA/time_series/$SAMPLE
    echo -e "bzip2 -9 $file" >> $log_file
    if [ -z $TESTING ]; then bzip2 -9 $file ; fi
    file=$file.bz2
    echo -e "mv $file $dir_arch/" >> $log_file
    if [ -z $TESTING ]; then mv $file $dir_arch/ ; fi
    echo >> $log_file
done

rmdir $dir_scratch

exit 0
