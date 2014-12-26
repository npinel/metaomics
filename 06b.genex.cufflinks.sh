#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 25 October 2013                                                                  # 
# Script: 06b.genex.cuffdiff.sh                                                          #
# Purpose:                                                                               #
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
  -a  gff annotation file (must correspond to reference sequence reads were aligned to)
  -o  output directory description (the run date is appended to the value provided)
  -b  subdirectory within the $LAOTS/sample/type base where the bam file is located (default NULL)
  -r  reference bam file if using Cuffdiff (not needed for Cufflinks)
  -i  description for bam file (e.g., if bam file is D49_MG.Bio17.x.clp7.mapq0.bam, use -i '_MG.Bio17.x.clp7.mapq0')
  -d  data type (MG = metagenomics; MT = metatranscriptomics)
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

####################
# define variables #
####################
FILE=
TYPE=MG
REFBAM=
BAMDESC=
OUTDIRDESC=
DATDIRBASE=
TESTING=0
SAVE_READS=0
while getopts "f:d:a:i:o:b:r:ht" OPTION
do
    case $OPTION in
	f) FILE=$OPTARG ;;
	d) TYPE=$OPTARG ;;
	a) ANNOT=$OPTARG ;;
	i) BAMDESC=$OPTARG ;;
	o) OUTDIRDESC=$OPTARG ;;
	b) DATDIRBASE=$OPTARG ;;
	r) REFBAM=$OPTARG ;;
	t) TESTING=1 ;;
	h) usage ;;
	?) usage ;;
    esac
done
DATDIRBASE='/'$DATDIRBASE
OUTDIRBASE='/cufflinks'
if [[ -n $OUTDIRDESC ]]; then
    OUTDIRBASE=$OUTDIRBASE.$OUTDIRDESC
fi
OUTDIRBASE=$OUTDIRBASE'.'`date +%Y%m%d`

# load ESB shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
export PATH=$PATH:$TOOLS

# module loading will only work within an oar job
# the same applies to the approach to determining number of cores available
# a number itself determined by the resources requested for the job
# modules last updated 26 October 2013
# check periodically for updates
# the current module version is not the latest Cufflinks version
#module load Cufflinks/2.0.2-goolf-1.4.10
#CUFFDIFF=`which cuffdiff`
#CUFFLINKS=`which cufflinks`

CUFFLINKS=$TOOLS/cufflinks-2.2.1.Linux_x86_64/cufflinks
CUFFDIFF=$TOOLS/cufflinks-2.2.1.Linux_x86_64/cuffdiff

if [ -z $CUFFLINKS ]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}An executable for Cufflinks cannot be found."
    exit 1
fi
if [ -z $CUFFDIFF ]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}An executable for Cuffdiff cannot be found."
    exit 1
fi
if [ -z $FILE ]; then
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldred}ERROR${bldblk}] ${txtrst}A sample list file, or a sample identifier are required."
    usage
fi

THREADS=`cat $OAR_RESOURCE_PROPERTIES_FILE | wc -l` # the file contains 1 line per core 
if [ -z $THREADS ]; then THREADS=1; fi

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

#######################
# execute the mapping #
#######################
for S in "${SAMPLES[@]}"
do
    DATDIR=$LAOTS/$S/$TYPE$DATDIRBASE
    OUTDIR=$LAOTS/$S/$TYPE$OUTDIRBASE
    BAM=$DATDIR/$S$BAMDESC.bam
    
    if [ -f $BAM ]; then
	if [ -z $REFBAM ]; then
	    echo -e "${bldblk}[`date +%H:%M:%S` ${bldgrn}START${bldblk}] ${txtrst}Starting to run Cufflinks on ${bldblk}$S${txtrst}."
	    $CUFFLINKS -p $THREADS -o $OUTDIR -G $ANNOT $BAM
	else
	    echo -e "${bldblk}[`date +%H:%M:%S` ${bldgrn}START${bldblk}] ${txtrst}Starting to run Cuffdiff on ${bldblk}$S${txtrst}."
	    $CUFFDIFF -u -p $THREADS -o $OUTDIR $ANNOT $REFBAM $BAM
	fi
	
	declare -a OUTFILES=($OUTDIR/*)
	for F in "${OUTFILES[@]}"
	do
	    if [[ -s $F ]]; then
		FNAME=`basename $F`
		FNAME=${S}.cufflinks$BAMDESC.$FNAME
		mv $F $OUTDIR/$FNAME
	    else
		rm $F
	    fi
	done
	echo -e "${bldblk}[`date +%H:%M:%S`   ${bldgrn}END${bldblk}] ${txtrst}Finished processing ${bldblk}$S${txtrst}."
    else
	echo -e "Cannot find the bam file $BAM"
    fi
done

exit 0
