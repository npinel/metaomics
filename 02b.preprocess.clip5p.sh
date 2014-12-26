#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 9 April 2013                                                                     #
# Script: 02b.preprocess.clip5p.sh                                                       #
# Purpose:                                                                               #
#    Performs a hard clip on the 5' of sequence reads to remove compositional biases     #
#    that may arise from instrument or sequencing adaptors.                              #
# NOTES:                                                                                 #
#    To establish the appropriate length to clip, examine compositional plots for DNA    #
#    sequence reads, not for reads derived from RNA libraries, as these appear inherent  #
#    ly uniform. Nevertheless, unless more information is obtained as to the effects     #
#    from RNA adaptors, it should be assumed that all reads suffer similar biases.       #
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
  -d  data type (MG=Metagenomic; MT=Metatranscriptomic; default MG)
  -p  position (base 1) for start of clipped read (e.g., -p 8 clips the first 7 bp)
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
POS=1
DATA_TYPE=MG
while getopts "f:p:d:ht" OPTION
do
    case $OPTION in
	h) usage ;;
	f) FILE=$OPTARG ;;
	p) POS=$OPTARG ;;
	d) DATA_TYPE=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done
let "CLP=$POS - 1"

if [[ -z $FILE ]]
then
    usage
fi

# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

# an intel-compiled version of fastx exists
# should probably add a function to search
# available modules and either interactively
# or automatically chose the desired/latest
# version - NP (15 Apr 2013)
#
# execution note: attempting through load the module this
# way through a passive oarsub submission fails
# consult how to invoke 'module load' in passive mode
# NP (15 Apr 2013)
#echo -e "loading fastx module: FASTX-Toolkit/0.0.13.2-ictce-4.0.6"
#`module load FASTX-Toolkit/0.0.13.2-ictce-4.0.6`

#######################
# identify input type #
#######################
SAMPLES=
if [ -f $FILE ]
then
    declare -a SAMPLES=()
    
    while read -r line; do
	[[ "$line" =~ ^#.*$ ]] && continue
	set -- $line
	declare -a fields=($*)
	SAMPLES=("${SAMPLES[@]}" ${fields[0]})
    done < $FILE
    
else
    declare -a SAMPLES=($FILE)
fi

############################
# loop through the samples #
############################
for S in "${SAMPLES[@]}"
do
    echo "Starting to process reads for $S"

    for R in 1 2 # may need to include an option to process single end as well
    do
	IFILE=$LAOTS/$S/$DATA_TYPE/${S}_$DATA_TYPE.R$R.fq
	if [ -e $IFILE ] ; then
	    OFILE=$LAOTS/$S/$DATA_TYPE/${S}_$DATA_TYPE.clp$CLP.R$R.fq
	    if [[ $TESTING -eq 1 ]]
	    then
		echo "testing: $TOOLS/fastx/fastx_trimmer -f $POS -Q 33 -i $IFILE -o $OFILE"
	    else
		echo "executing: $TOOLS/fastx/fastx_trimmer -f $POS -Q 33 -i $IFILE -o $OFILE"
		`$TOOLS/fastx/fastx_trimmer -f $POS -Q 33 -i $IFILE -o $OFILE`
	    fi
	fi
    done
done

exit 0
