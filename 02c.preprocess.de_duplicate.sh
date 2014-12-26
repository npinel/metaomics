#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 8 April 2013                                                                     #
# Script: 02c.preprocess.de_duplicate.sh                                                 #
# Purpose:                                                                               #
#    To remove suspected duplicate reads from metagenomic files using fastuniq           #
#    (http://sourceforge.net/projects/fastuniq/).                                        #
#    Unless a good known reason exists for doing so, duplicates shouldn't be removed     #
#    from metatranscriptomic reads, since bona fide exact copies of pairs,               #
#    given the nature of transcripts, are likely to be found.                            #
#                                                                                        #
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
  -d  data type; values = MG|MT (required)
  -f  file with sample list, or single dataset name (required)
  -s  suffix characterizing the read file name (e.g., clp7)
  -t  test, flag
  -h  display this message

EOF
    exit 1
}

################
# LOAD MODULES #
################

## no modules available

####################
# DEFINE VARIABLES #
###################
# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

WD=`pwd`
FILE=
TYPE=
SUFF=
TESTING=

while getopts "d:f:s:th" OPTION
do
    case $OPTION in
	h) usage ;;
	d) TYPE=$OPTARG ;;
	f) FILE=$OPTARG ;;
	s) SUFF=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

if [[ -z $FILE ]] || [[ -z $TYPE ]]
then
    usage
fi

############################
# identify input file type #
############################

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

for S in ${SAMPLES[@]} ; do
    R1=$LAOTS/$S/$TYPE/${S}_$TYPE.$SUFF.R1.fq
    R2=$LAOTS/$S/$TYPE/${S}_$TYPE.$SUFF.R2.fq
    echo "Checking for the existence of $R1 and $R2"
    if [ -e $R1 ] && [ -e $R2 ] ; then
	# create the list file
	echo -e "$R1\n$R2" > "$LAOTS/$S/$TYPE/fastuniq.list"
	R1=$LAOTS/$S/$TYPE/${S}_$TYPE.$SUFF.uniq.R1.fq
	R2=$LAOTS/$S/$TYPE/${S}_$TYPE.$SUFF.uniq.R2.fq
	# execute the command
	echo "Executing: $TOOLS/fastuniq -i $LAOTS/$S/$TYPE/fastuniq.list -o $R1 -p $R2"
	`$TOOLS/fastuniq -i $LAOTS/$S/$TYPE/fastuniq.list -o $R1 -p $R2`
	rm $LAOTS/$S/$TYPE/fastuniq.list
	echo "Finished: $S"
    fi
done

exit 0
