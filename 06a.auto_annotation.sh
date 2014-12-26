#!/bin/bash -l
##########################################################################################
# Author(s):                                                                             #
#    1. Shaman Narayanasamy (shaman.narayanasamy@uni.lu)                                 #
# Affiliation(s):                                                                        #
#    1. Luxembourg Center for Systems Biomedicine/University of Luxembourg               #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 17 May 2013                                                                      #
# Script: 06a.auto_annotation.sh                                                         #
# Purpose:                                                                               #
#         Runs gene calling/prediction softwares on assemblies and raw sequencing reads. #
#	  Output will be annotations. Output formats will be determined by the type of   #
#	  program used.									 #
# NOTES:                                                                                 #
#    As of Apr 17, it is still in progress. Using Nic's shell scripts as templates.      #
#      1. README files available for each of the indivudual gene prediction softwares.   # 
#      2. FragGeneScan is not working.                                                   #  
#      3. Prodigal is suitable to be used on isolate genomes, while MetaGeneMark may be  #
#         more suitable for the metagenomic assemblies.                                  #
##########################################################################################
## INITIALIZE FUNCTIONS
##########################################################################################
## function: usage
## purpose: if you don't know how the script works

usage(){
    cat <<EOF

Usage: $0 options

OPTIONS:
  -t  testing = 1|0
  -p  program; values = pdg (prodigal)|mgm (MetaGeneMark)|fgs (FragGeneScan) | all 
  -i  input file
  -h  display this message

EOF
exit 1
}

##########################################################################################
## function: run prodigal
## purpose runs prodigal... DUH!
## input: fasta file..?
## output: gbk (default), gff, or sco. Nucleotide and/or protein fasta (optional)

run_pdg(){
  testing=$1
  INPUT=$2
  OUTPUT=`echo "$INPUT" | sed 's/fasta/gff/'` # Output file name
  FAA=`echo "$INPUT" | sed 's/fasta/faa/'` # Amino acid fasta output
  FNA=`echo "$INPUT" | sed 's/fasta/fna/'` # Nucleotide fasta output
  if [ $testing -eq 0 ]
  then
    echo "Executing: prodigal -a -d -f gff -i $INPUT -o $OUTPUT"
    echo "Producing: $OUTPUT, $FAA and $FNA"
    `prodigal -a $FAA -d $FNA -f gff -i $INPUT -o $OUTPUT`
    # -a: Writes protein sequence gene file
    # -d: Writes nucleotide gene file
    # -f: Annotation file format
    # -o: Out put file format
    # NOTE: For other options, run prodigal and modify this function accordingly
  elif [ $testing -eq 1 ]
  then
    echo "Testing: prodigal -a -d -f gff -i $INPUT -o $OUTPUT"
    echo "Producing: $OUTPUT, $FAA and $FNA"
  else [ -z $INPUT ]
    usage
  fi
}

##########################################################################################
## function: run MetaGeneMark
## purpose calls genes with MetaGeneMark
## input: fasta file
## output: list, gff, nucleotide and/or protien fasta (optional)
## NOTE: More commands can be added to generate fasta files for the protein/nucleotide
##       sequences.

run_mgm(){
  testing=$1
  INPUT=$2
  OUTPUT=`echo "$INPUT.gff"`
  MODEL=/mnt/nfs/projects/ecosystem_biology/local_tools/genePredictionTools/MetaGeneMark_linux64/MetaGeneMark_v1.mod
  if [ $testing -eq 0 ]
  then
    echo "Executing: gmhmmp -a -d -f G -m MetaGeneMark_v1.mod -o $OUTPUT $INPUT"
    `gmhmmp -a -d -f G -m $MODEL -o $OUTPUT $INPUT`
    # -o: Output file name
    # -a: Writes protein sequence gene file
    # -d: Writes nucleotide sequence gene file
    # -f: Output format (set as GFF2)
  elif [ $testing -eq 1 ]
  then
    echo "Testing: gmhmmp -a -d -f G -m $MODEL -o $OUTPUT $INPUT"
  else [ -z $INPUT ]
    usage
  fi
}

##########################################################################################
## function: run FragGeneScan
## purpose runs FragGeneScan... DUH!
## input: fasta file or sequencing reads
## output: gene coordianate list, nucleotide fasta and protein fasta
## NOTE: Was not working previously, tried again and it seems to work well

run_fgs(){
  testing=$1
  INPUT=$2
  OUTPUT=`echo "$INPUT.gff"`
  if [ $testing -eq 0 ]
  then
    echo "executing: fgs -genome=$INPUT -out=$OUTPUT -complete=1 -train=complete"
    `fgs -genome=$INPUT -out=$OUTPUT -complete=1 -train=complete`
    # -genome: genome consensus sequence fasta file
    # -out: output file name
    # -complete: 1=complete genome seq|0=sequencing reads
    # -train: training models (in train directory of FGS)
  elif [ $testing -eq 1 ]
  then
    echo "testing: fgs -genome=$INPUT -out=$OUTPUT -complete=1 -train=complete"
  else [ -z $INPUT ]
    usage
  fi
}

##########################################################################################
## function: run CD-HIT
## purpose: Merges multiple faa (amino acid) sequence files
## input: list of faa files
## output: merged faa file that can be used for annotation. later on
## NOTE: To be written

##########################################################################################
## function: merge gff files
## purpose: Merges multiple gff files
## input: list of gff files to merge
## output: merged faa file that can be used for annotation. later on
## NOTE: To be written


####################
# DEFINE VARIABLES #
####################

. /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
shopt -s expand_aliases # Enables aliases. Since the programs are aliases in the
                        # esb.shortcuts file.
INDIR=`pwd`
PROG=
INPUT=
TESTING=

while getopts "ht:p:i:" OPTION
do
  case $OPTION in
    h) usage ;;
    p) PROG=$OPTARG ;;
    i) INPUT=$OPTARG ;;
    t) TESTING=$OPTARG ;;
    ?) usage ;;
  esac
done
####################

#################
# RUN FUNCTIONS #
#################

if [[ "$PROG" = "pdg" ]]
then
  echo "Gene prediction using Prodigal"
  run_pdg $TESTING $INPUT
elif [[ "$PROG" = "mgm" ]]
then
  echo "Gene prediction using MetaGeneMark"
  run_mgm $TESTING $INPUT
elif [[ "$PROG" = "fgs" ]]
then
  echo "Gene prediction using FragGeneScan"
  run_fgs $TESTING $INPUT
elif [[ "$PROG" = "all" ]]
then
  echo "Gene prediction using Prodigal, FragGeneScan and Metagenemark"
  run_fgs $TESTING $INPUT
  run_mgm $TESTING $INPUT
  run_pdg $TESTING $INPUT
else [ -z $PROG ]
  echo "Error: Gene prediction program not specified"
  usage
fi

#################

exit 0
