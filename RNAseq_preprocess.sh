#!/bin/bash -l
####################################################################################
# Author(s):                                                                       #
#    1. Shaman Narayanasamy (shaman.narayanasamy@uni.lu)                           #
# Affiliation(s):                                                                  #
#    1. Luxembourg Center for Systems Biomedicine/University of Luxembourg         #
# Project(s):                                                                      #
#    Any RNA-seq sequencing data sets					           #
# Date: 15 July 2013                                                               #
# Script: RNAseq_preprocess.sh							   #
# Purpose: Preprocessing pipeline for RNAseq reads                                 #
# NOTES: As of 15th July                                                           #
#      Script runs AdapterRemoval program on a given paired end RNAseq data set.   #
#      Other functions will be added soon to include quality checks and other	   #
#      quality filtering programs for more robust preprocessing                    #
####################################################################################
## INITIALIZE FUNCTIONS
####################################################################################
## function: usage
## purpose: if you don't know how the script works

usage(){
    cat <<EOF

Usage: $0 options
    
OPTIONS:
    -t testing = 1|0
    -1 fastq1
    -2 fastq2
    -o output name (only basename)
    -h print this help message

EOF
exit 1
}

####################################################################################
## function: runs AdapterRemoval on paired end data
## purpose: adapter trimming, quality filtering and **lateral collapsing
## input: a pair of FASTQ files (uncompressed)
## output: 1) **collapsed; 2) single; 3) pairs1; 4) pairs2; 5) discarded

run_adapterRemoval(){
  testing=$1
  FQ1=$2
  FQ2=$3
  BASE=$4

  if [ $testing -eq 1 ]
  then
    echo "testing: AdapterRemoval --file1 $FQ1 --file2 $FQ2 --trimns --trimqualities --minquality 30 --5prime AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --basename $BASE --collapse"
  elif [ $testing -eq 0 ]
  then
    echo "Starting for pairs $FQ1 and $FQ2 on `date`"

    echo "executing: AdapterRemoval --file1 $FQ1 --file2 $FQ2 --trimns --trimqualities --minquality 30 --5prime AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --basename $BASE --collapse"
    /home/users/snarayanasamy/AdapterRemoval-1.5/AdapterRemoval --file1 $FQ1\
      --file2 $FQ2\
      --trimns --trimqualities --minquality 30\
      --5prime AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG\
      --basename $BASE\
      --collapse

    echo "Completed for pairs $FQ1 and $FQ2 on `date`"
  else [ -z $FQ1 ]
    usage
  fi
}

####################################################################################
## function: decompresses gunzipped files into given directory
## purpose: decompress gunzipped files
## input: a pair of FASTQ files
## output: unzipped files

run_gunzip(){
	testing=$1
	FQ1=$2
	FQ2=$3
	OUTDIR=$4
	OUTPUT1=`echo $FQ1 | rev | cut -f1 -d "/" | rev | sed 's/.gz//'`
	OUTPUT2=`echo $FQ2 | rev | cut -f1 -d "/" | rev | sed 's/.gz//'`
	if [ $testing -eq 1 ]
	then
		echo "testing: Unzipping files..."
		echo "testing: gunzip -c $FQ1 > $OUTDIR/$OUTPUT1"
		echo "testing: gunzip -c $FQ2 > $OUTDIR/$OUTPUT2"
	elif [ $testing -eq 0 ]
	then
		echo "Unzipping files..."
		echo "executing: gunzip -c $FQ1 > $OUTDIR/$OUTPUT1"
		gunzip -c $FQ1 > $OUTDIR/$OUTPUT1
		echo "executing: gunzip -c $FQ2 > $OUTDIR/$OUTPUT2"
		gunzip -c $FQ2 > $OUTDIR/$OUTPUT2
	else [ -z $FQ1 ]
		usage
	fi
}

####################################################################################
## function: decompresses and untars fastq files into given directory
## purpose: decompress and untar fastq files
## input: a pair of FASTQ files
## output: unzipped files

run_untar(){
	testing=$1
	FQ1=$2
	FQ2=$3
	OUTDIR=$4
	if [ $testing -eq 1 ]
	then
		echo "testing: Untarring files..."
		echo "testing: tar -C $OUTDIR -zxvf $FQ1"
		echo "testing: tar -C $OUTDIR -zxvf $FQ2"
	elif [ $testing -eq 0 ]
	then
		echo "Untarring files..."
		echo "executing: tar -C $OUTDIR -zxvf $FQ1"
		tar -C $OUTDIR -zxvf $FQ1
		echo "executing: tar -C $OUTDIR -zxvf $FQ1"
		tar -C $OUTDIR -zxvf $FQ2
	else [ -z $FQ1 ]
		usage
	fi
}

####################################################################################

OUTDIR=`pwd`
FQ1=
FQ2=
OUTPUT=
TESTING=0

while getopts "ht:1:2:o:" OPTION
do
  case $OPTION in
    h) usage ;;
    1) FQ1=$OPTARG ;;
    2) FQ2=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    t) TESTING=$OPTARG ;;
    ?) usage ;;
  esac
done

if [[ $FQ1 =~ ".tar.gz" ]]
then
  run_untar $TESTING $FQ1 $FQ2 $OUTDIR
  FQ1_2=`echo $FQ1 | rev | cut -f1 -d "/" | rev | sed 's/.tar.gz//'`
  FQ2_2=`echo $FQ2 | rev | cut -f1 -d "/" | rev | sed 's/.tar.gz//'`
  run_adapterRemoval $TESTING $OUTDIR/$FQ1_2 $OUTDIR/$FQ2_2 $OUTPUT
elif [[ $FQ1 =~ ".gz" ]]
then
  run_gunzip $TESTING $FQ1 $FQ2 $OUTDIR
  FQ1_2=`echo $FQ1 | rev | cut -f1 -d "/" | rev | sed 's/.gz//'`
  FQ2_2=`echo $FQ2 | rev | cut -f1 -d "/" | rev | sed 's/.gz//'`
  run_adapterRemoval $TESTING $OUTDIR/$FQ1_2 $OUTDIR/$FQ2_2 $OUTPUT
else
  run_adapterRemoval $TESTING $FQ1 $FQ2 $OUTDIR/$OUTPUT
fi

exit 0
