#!/bin/bash
##########################################################################################
# Author(s):                                                                             #
#    1. Nicolas Pinel (npinel@systemsbiology.org)                                        #
# Affiliation(s):                                                                        #
#    1. Institute for Systems Biology / Luxembourg Center for Systems Biomedicine        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 10 April 2014                                                                    # 
# Script: 04c.analysis.metaphlan.sh                                                      #
# Purpose:                                                                               #
#    Bash shell wrapper script to conduct MetaPhlAn analyses                             #
# NOTES:                                                                                 #
#    Full description of MetaPhlAn (Segata et al., 2012), is available from:             #
#    https://bitbucket.org/nsegata/metaphlan                                             #
##########################################################################################
# oarsub directives
#OAR --stdout log.oarjob_%jobid%.bwa.txt
#OAR --stderr log.oarjob_%jobid%.bwa.txt
############################
#
# define message colors
# usage:
# echo -e "${bldred}<stylized text here>${txtrst}<non-stylized text here>"
bldblk='\e[1;30m' # bold black
bldred='\e[1;31m' # bold red
bldgrn='\e[1;32m' # bold green
bldblu='\e[1;34m' # bold blue
txtbld='\e[1m'    # bold, default color
txtrst='\e[0m'    # reset text style

# load ESB shortcuts
if [[ -z $ESB ]]; then
    source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
    export PATH=$PATH:$TOOLS
fi

# load needed modules
module load Python/2.7.6-goolf-1.4.10 # needed if interleaving files
module load Bowtie2/2.0.2-ictce-5.3.0
module load BLAST/2.2.28-ictce-5.3.0-Python-2.7.3

# define functions
usage(){
    cat << EOF

Usage: $0 options

OPTIONS:
  -f  file with tab-separated sample list (ID\tYYYY-MM-DD; required)
  -d  data type (MG vs MT; default=MG)
  -q  qualifier string for read files (e.g., if read file is D11_MG.clp7.no_ovlp.qv20.lng40.R1.fq,
      the qualifier corresponds to 'clp7.uniq.qv30.ln40')
  -t  testing flag; prints out the commands that would have been executed
  -h  display this message

EOF
    exit 1
}

####################
# define variables #
####################
FILE=
DATA=MG
QUAL=
TESTING=0

while getopts "f:q:d:ht" OPTION
do
    case $OPTION in
	h) usage ;;
	f) FILE=$OPTARG ;;
	q) QUAL=$OPTARG ;;
	d) DATA=$OPTARG ;;
	t) TESTING=1 ;;
	?) usage ;;
    esac
done

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
    echo -e "${bldblk}[`date +%H:%M:%S`  ${bldgrn}START${bldblk}] ${txtrst}Processing ${bldblk}$S${txtrst}."
    DATDIR=$LAOTS/$S/$DATA
    R1=$DATDIR/${S}_$DATA.$QUAL.R1.fq
    R2=$DATDIR/${S}_$DATA.$QUAL.R2.fq
    RO=$DATDIR/${S}_$DATA.$QUAL.R12.fq
    BO=$DATDIR/$S.metaphlan.$QUAL.R12.bt2out.txt # bowtie2 output
    CO=$DATDIR/$S.metaphlan.$QUAL.R12.classification.txt # classification output
    # interleave sequence reads
    echo -e "${bldblk}[`date +%H:%M:%S`   ${bldblu}PROG${bldblk}] ${txtrst}Interleaving sequences."
    /home/users/npinel/sandbox/interleave_seqs.py -r1 $R1 -r2 $R2 -o $RO
    # run MetaPhlAn
    echo -e "${bldblk}[`date +%H:%M:%S`   ${bldblu}PROG${bldblk}] ${txtrst}Running MetaPhlAn."
    $TOOLS/metaphlan.py $RO --bowtie2db $TOOLS/metaphlan/bowtie2db/mpa --bt2_ps sensitive-local --nproc $THREADS --bowtie2out $BO > $CO
    # clean up
    echo -e "${bldblk}[`date +%H:%M:%S`   ${bldblu}PROG${bldblk}] ${txtrst}Cleaning up."
    rm $RO
    echo -e "${bldblk}[`date +%H:%M:%S` ${bldgrn}FINISH${bldblk}] ${txtrst}Finished processing ${bldblk}$S${txtrst}."
done

exit 0
