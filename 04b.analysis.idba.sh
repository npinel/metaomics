#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 25 October 2013                                                                  # 
# Script: 04b.analysis.idba.sh                                                           #
# Purpose:                                                                               #
#    Bash shell wrapper script to configure and run idba assemblies                      #
# NOTES:                                                                                 #
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

usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -d  data type; values = Metagenomic|MG|Metatranscriptomic|MT (required)
  -f  file with sample list, or single dataset name (required)
  -o  output directory suffix
  -q  qualifier string
  -h  display this message

EOF
    exit 1
}

####################
# DEFINE VARIABLES #
####################
FILE=
DATA=MG
QUAL=
ODIRSUFF=
CLEAN=
TESTING=0
MAPQ=5
while getopts "f:d:o:q:htc" OPTION
do
    case $OPTION in
	f) FILE=$OPTARG ;;
	d) DATA=$OPTARG ;;
	o) ODIRSUFF=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	t) TESTING=1 ;;
	c) CLEAN=1 ;;
	h) usage ;;
	?) usage ;;
    esac
done

# load ESB shortcuts
if [[ -z $ESB ]]; then
    source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
    export PATH=$PATH:$TOOLS
fi

if [[ -z $FILE ]]; then
    echo -e "${bldblk}[${bldred}ERROR\t${bldblk}`date +%H:%M:%S`]  ${txtrst}A sample list file, or a sample identifier are required."
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

################################
# define additional parameters #
################################
BDIR=$TOOLS/idba-1.1.1.icc/bin # as of 2013-10-31, no modules exist for idba
mink=25
maxk=95
step=5
PERID=0.98 # percent sequence identity to compbine contigs
THREADS=`cat $OAR_RESOURCE_PROPERTIES_FILE | wc -l` # the file contains 1 line per core
if [ -z $THREADS ]; then THREADS=1; fi

for S in "${SAMPLES[@]}"
do
    DDIR=${LAOTS}/${S}/$DATA
    WDIR=${LAOTS}/${S}/$DATA
    ODIR=$WDIR/$ODIRSUFF
    FA=$WDIR/${S}_$DATA.${QUAL}.PAIRS.fa
    # convert fastq to fasta and merge if necessary
    if [ ! -f $FA ] ; then
	echo -e "${bldblk}[PROG\t`date +%H:%M:%S`]  ${txtrst}Preparing reads data file for $S."
	echo -e "> $BDIR/fq2fa --merge $WDIR/${S}_$DATA.$QUAL.R1.fq $WDIR/${S}_$DATA.$QUAL.R2.fq $FA"
    	`$BDIR/fq2fa --merge $WDIR/${S}_$DATA.$QUAL.R1.fq $WDIR/${S}_$DATA.$QUAL.R2.fq $FA`
    fi
# Runs idba_tran for metatranscriptomic data sets. User
# specifies "MT" as a parameter
#    if [[ $DATA -eq 'MT' ]]; then
    if [[ $DATA == 'MT' ]]; then
	echo -e "${bldblk}[${bldgrn}START${bldblk}\t`date +%H:%M:%S`]  ${txtrst}Starting assembly of metatranscriptomic data for $S."
	echo -e "> $BDIR/idba_tran -r $FA -o $ODIR --mink $mink --maxk $maxk --step $step --num_threads $THREADS --similar 0.98 --pre_correction"
    	`$BDIR/idba_tran -r $FA -o $ODIR --mink $mink --maxk $maxk --step $step --num_threads $THREADS --similar $PERID --pre_correction`

# Runs idba_ud for metagenomic data sets. User specifies
# "MG" as a parameter
### BUG: The script does not seem to be able to run idba_ud.
###	 It automatically does idba_tran
    else
	echo -e "${bldblk}[${bldgrn}START${bldblk}\t`date +%H:%M:%S`]  ${txtrst}Starting assembly of metagenomic data for $S."
	echo -e "> $BDIR/idba_ud -r $FA -o $ODIR --mink $mink --maxk $maxk --step $step --num_threads $THREADS --pre_correction"
    	`$BDIR/idba_ud -r $FA -o $ODIR --mink $mink --maxk $maxk --step $step --num_threads $THREADS --similar $PERID --pre_correction`
    fi

    mv $ODIR/contig.fa $WDIR/${S}_${DATA}.$ODIRSUFF.fa
    mv $ODIR/log $WDIR/${S}_${DATA}.$ODIRSUFF.log
    if [[ $CLEAN -eq 1 ]]; then
	rm -fr $ODIR
    fi
    echo -e "${bldblk}[${bldgrn}END${bldblk}\t`date +%H:%M:%S`]  ${txtrst}Finished assembling $S."
done

exit 0
