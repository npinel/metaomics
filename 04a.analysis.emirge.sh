#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 22 August 2013                                                                   #
# Script: 04a.analysis.emirge.sh                                                         #
# Purpose:                                                                               #
#                                                                                        #
##########################################################################################

# define message colors
# usage: echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
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
  -s  suffix characterizing the input/output files (required)
  -t  test, flag
  -h  display this message

EOF
    exit 1
}

####################
# DEFINE VARIABLES #
###################
# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
source /home/users/npinel/.profile

FILE=
TYPE=
SUFF=
TESTING=
DB_FILE=/mnt/nfs/projects/ecosystem_biology/databases/rRNA/SSU/SILVA_SSU_111_NR_tax_silva_trunc_dna.fasta
DB_BASE=/mnt/nfs/projects/ecosystem_biology/databases/rRNA/SSU/SILVA_SSU_111_NR_tax_silva_trunc_dna
ITER=100
THREADS=12
INSERT_AVE=400
INSERT_SD=100
MAX_LENGTH=100

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
ODIR=emirge.$SUFF

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
	WDIR=$LAOTS/$S/$TYPE/$ODIR
	CMD="$TOOLS/emirge.py $WDIR -1 $R1 -2 $R2 -f $DB_FILE -b $DB_BASE -l $MAX_LENGTH -i $INSERT_AVE -s $INSERT_SD -n $ITER -a $THREADS --phred33"
	echo "Executing: $CMD"
	$CMD

	## now process output
	for i in 50 60 70 80 90 100
	do
	    dir="$WDIR/iter.$i"
	    out="$WDIR/${S}_$TYPE.$SUFF.emirge_iter$i.fasta"
	    echo -e "Executing: $TOOLS/emirge_rename_fasta.py $dir > $out"
	    $TOOLS/emirge_rename_fasta.py $dir > $out
	done

    fi
done

exit 0
