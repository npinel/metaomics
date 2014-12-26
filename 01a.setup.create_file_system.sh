#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 9 April 2013                                                                     #
# Script: 01a.setup.create_file_system.sh                                                #
# Purpose:                                                                               #
#    Create archive and work directory structure for new samples from the 1-year         #
#    time series LAO project.                                                            #
# NOTES:                                                                                 #
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
  -f  file with tab-separated sample list (ID\tYYYY-MM-DD; required)
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

####################
# DEFINE VARIABLES #
####################
PWD=`pwd`
FILE=
TESTING=
while getopts "f:ht" OPTION
do
    case $OPTION in
	h) usage ;;
	f) FILE=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

if [[ -z $FILE ]]
then
    usage
fi

source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

################################################
# reads sample list and loop through the files #
################################################
while read -r line; do
    [[ "$line" =~ ^#.*$ ]] && continue

    set -- $line
    declare -a fields=($*)
    sample_id=${fields[0]}
    sample_date=${fields[1]}
    a_mg=$ARCHIVE/LAO/Metagenomic
    a_mt=$ARCHIVE/LAO/Metatranscriptomic
    work_dir=$LAO/$sample_date
    w_mg=$work_dir/Metagenomic
    w_mt=$work_dir/Metatranscriptomic

    if [[ $TESTING -eq 1 ]]; then
	if [[ ! -d "$work_dir" ]]; then
	    echo -e "\nCreating directory structure for $sample_id"
        # create directories in $ARCHIVE
	    echo "cd $a_mg; mkdir $sample_date; cd one_year_time_series; ln -s ../$sample_date $sample_id"
	    echo "cd $a_mt; mkdir $sample_date; cd one_year_time_series; ln -s ../$sample_date $sample_id"
        # create directories in $LAO
	    echo "mkdir $work_dir"
	    echo "mkdir $w_mg; mkdir $w_mt; cd $work_dir; ln -s Metagenomic MG; ln -s Metatranscriptomic MT"
        # create symlinks in $LAOTS
	    echo "cd $LAOTS; ln -s ../$sample_date $sample_id"
        # return to invocation directory
	    echo "cd $PWD"
	else
	    echo -e "\nDirectories for $sample_id already exist. OMITTING."
	fi
    else
	if [[ ! -d "$work_dir" ]]; then
	    `cd $a_mg; mkdir $sample_date; cd one_year_time_series; ln -s ../$sample_date $sample_id`
	    `cd $a_mt; mkdir $sample_date; cd one_year_time_series; ln -s ../$sample_date $sample_id`
	    `mkdir $work_dir`
	    `mkdir $w_mg; mkdir $w_mt; cd $work_dir; ln -s Metagenomic MG; ln -s Metatranscriptomic MT`
	    `cd $LAOTS; ln -s ../$sample_date $sample_id`
	    `cd $PWD`
	fi
    fi

done < $FILE

exit 0
