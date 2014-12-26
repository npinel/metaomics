#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 20 February 2014                                                                 # 
# Script: util.batch_command.sh                                                          #
# Purpose:                                                                               #
#    Input a command to be applied to all samples                                        #
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
  -c  command
  -t  test command without executing
  -h  display this message

NOTES:
  1) This script will change directory to $LAOTS prior to executing, thus, the environmental
variable $LAOTS should not be included in the command string. Upon execution of all the commands,
the script will return to the invocation directory.
  2) To indicate the locations in the command string that should be replaced by the sample name,
use %s. E.g.,

~/git/laoPipeline/util.batch_command.sh [-t] -f sample_list.txt -c 'mv %s/MT/%s_MT.old_file.fq %s/MT/%s_MT.new_file.fq'

EOF
    exit 1
}

####################
# define variables #
####################
FILE=
COMMAND=
TESTING=
while getopts "f:c:ht" OPTION
do
    case $OPTION in
	f) FILE=$OPTARG ;;
        c) COMMAND=$OPTARG ;;
	t) TESTING=1 ;;
	h) usage ;;
	?) usage ;;
    esac
done

if [ -z $FILE ]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}A sample list file, or a sample identifier are required."
    usage
fi

source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

#################################################
# identify input type and populate sample array #
#################################################
declare -a SAMPLES=()
while read -r line; do
    [[ "$line" =~ ^#.*$ ]] && continue
    set -- $line
    declare -a fields=($*)
    SAMPLES=("${SAMPLES[@]}" ${fields[0]})
done < $FILE

NUM=`echo $COMMAND | grep -o -e '%' | wc -l`

PWD=`pwd`
cd $LAOTS
###########
# execute #
###########
for S in "${SAMPLES[@]}"
do
    ARG=
    for (( c=1; c<=$NUM; c++ )); do ARG=$ARG" $S"; done
    printf -v CMD "$COMMAND"$ARG

    if [[ $TESTING -eq 1 ]]; then
	echo $CMD
    else
	echo $CMD
	$CMD
    fi
done
cd $PWD

exit 0
