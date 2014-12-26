#!/bin/bash -l
##########################################################################################
# Author(s):                                                                             #
#    1. Shaman Narayanasamy (shaman.narayanasamy@uni.lu)                                 #
#    2. Nic Pinel (npinel@systemsbiology.org)                                            #
# Affiliation(s):                                                                        #
#    1. Luxembourg Center for Systems Biomedicine/University of Luxembourg               #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 17 April 2013                                                                    #
# Script: 05c.variant_identification.sh                                                  #
# Purpose:                                                                               #
#    Calls variants for the Microthrix and Zooglea alignements. Outputs variant call     #
#    files (.vcf) and consensus fastq file.                                              # 
# NOTES:                                                                                 #
#    As of Apr 17, it is still in progress. Using Nic's shell scripts as templates.      #
#    Runs three variant calling programs; i) samtools mpileup				 #
#				          ii) freebayes					 #
#					  iii) GATK UnifiedGenotyper **			 #
#    ** Special preprocessing of bam files required before running the variant calls     #
#    Might include a final step of merging all the outputs from the different variant    #
#    callers.										 #
#											 #
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
  -i  input file
  -o  output type; values = vcf(default) |fastq|fasta
  -r  reference genome/metagenome = cm/cm_old (Candidatus microthrix)|zg (Zooglea)|mg (Metagenome)
  -h  display this message

EOF
exit 1
}

##########################################################################################
## function: generates VCF file using the mpileup pipeline
## purpose: runs samtools mpileup and pipes the output to bcftools and vcfutils
##          to generate a final vcf output
## input: SORTED BAM format file (binary)
## output: BCF format file (binary)
## NOTE: I merged the vcfutils function into this one, so that mpileup runs as a
##       single pipeline

run_mpileup(){
	testing=$1
	REF=$2
	BAM=$3
	VCF=`echo $BAM | rev | cut -f1 -d "/" | rev | sed 's/.bam/.mpu.vcf/'`
	# Cuts off the path and renames the file with the output extention. Maybe
	# there is a better way to do it...

	if [ $testing -eq 1 ]
	then
		echo "testing: samtools mpileup -uf $REF $BAM | bcftools view -vcg - | vcfutils.pl varFilter -d10 -D100 > $OUTDIR/$VCF"
	elif [ $testing -eq 0 ]
	then
		echo "executing: samtools mpileup -uf $REF $BAM | bcftools view -vcg - | vcfutils.pl varFilter -d10 -D100 > $OUTDIR/$VCF"
		samtools mpileup -uf $REF $BAM | bcftools view -vcg - | vcfutils.pl varFilter -d10 -D100 > $OUTDIR/$VCF
      else [ -z $BAM ]
	usage
      fi
}

##########################################################################################
## function: generates vcf file using freebayes variant caller
## purpose: run freebayes
## input: sorted bam file
## output: vcf file

run_freebayes(){
	testing=$1
	REF=$2
	BAM=$3
	VCF=`echo $BAM | rev | cut -f1 -d "/" | rev |sed 's/.bam/.frb.vcf/'`

	if [ $testing -eq 1 ]
	then 
		echo "testing freebayes -f $REF $BAM > $OUTDIR/$VCF"
	elif [ $testing -eq 0 ]
	then
		echo "executing freebayes -f $REF $BAM > $OUTDIR/$VCF"
		/mnt/nfs/projects/ecosystem_biology/local_tools/gkno_launcher/tools/freebayes/bin/freebayes -f $REF $BAM > $OUTDIR/$VCF
	else [ -z $BAM ]
		usage
	fi
}

##########################################################################################
## function: generates vcf file using GATK's UnifiedGenotyper variant caller
## purpose: runs UnifiedGenotyper
## input: sorted bam file with corresponding indexed bam within same folder
## output: vcf file
## NOTE: We should consider fitting in all the GATK analysis pipeline steps within this
##       Script

run_gatk(){
	testing=$1
	REF=$2
	BAM=$3
	VCF=`echo $BAM | rev | cut -f1 -d "/" | rev |sed 's/.bam/.gatk.vcf/'`
	if [ $testing -eq 1 ]
	then 
		echo "testing java -Xmx6g -jar /mnt/nfs/projects/ecosystem_biology/local_tools/gatk.jar -T UnifiedGenotyper -I $BAM -R $REF -o $OUTDIR/$VCF -glm BOTH"

	elif [ $testing -eq 0 ]
	then
		echo "executing java -Xmx6g -jar /mnt/nfs/projects/ecosystem_biology/local_tools/gatk.jar -T UnifiedGenotyper -I $BAM -R $REF -o $OUTDIR/$VCF -glm BOTH"
		java -Xmx6g -jar /mnt/nfs/projects/ecosystem_biology/local_tools/gatk.jar -T UnifiedGenotyper -I $BAM -R $REF -o $OUTDIR/$VCF -glm BOTH
	else [ -z $BAM ]
		usage
	fi
}

##########################################################################################
## function: generates consensus fastq file from from bam file
## purpose: if consensus fastq file is required
## input: BCF file (binary)
## output: consensus FASTQ file (binary)

generate_fastq(){
	testing=$1
	REF=$2
	BAM=$3
	FASTQ=`echo $BAM | rev | cut -f1 -d "/" | rev | sed 's/.bam/.fq/'`
	#FASTQ=`echo "$BAM" | sed 's/bam/fq/'` # Output file name
		
	if [ $testing -eq 1 ]
	then
		echo "testing: samtools mpileup -uf $REF $BAM | bcftools view -cg - | vcfutils.pl vcf2fq > $OUTDIR/$FASTQ"
	elif [ $testing -eq 0 ]
	then
		echo "executing samtools mpileup -uf $REF $BAM | bcftools view -cg - | vcfutils.pl vcf2fq > $OUTDIR/$FASTQ"
		samtools mpileup -uf $REF $BAM | bcftools view -cg - | vcfutils.pl vcf2fq > $OUTDIR/$FASTQ
	else [ -z $BAM ]
	  usage
	fi
}

##########################################################################################
## FUNCTIONS TO BE ADDED LATER ON
##########################################################################################

##########################################################################################
## function: generates consensus fastq file from from bam file
## purpose: if consensus fastq file is required
## input: consensus FASTQ file (text)
## output: consensus FA or FASTA file
## NOTE: Having a hard time finding a good fastq to fasta converter. Bloody hell!
## Please let me know if you have any reliable ones. Functions to generate fasta
## files will be unavailable until we solve this problem...
#
#generate_fasta(){
#	testing=$1
#	FASTQ=$2
#	FASTA=`echo "$FASTQ" | sed 's/fastq/fa/'` # Output file name
#	if [ $testing -eq 1 ]
#	then
#	  echo "testing: perl $TOOLS/fastq2fasta.pl $FASTQ"
#	elif [ $testing -eq 0]
#	then 
#	  echo "executing: perl $TOOLS/fastq2fasta.pl $FASTQ"
#	  perl $TOOLS/fastq2fasta.pl 
#	else [ -z $FASTQ]
#	  usage
#	fi
#}
#
##########################################################################################
## function: reports variant statistics for both SNPs and indels
## purpose: produce simple stats and plots for population variation
## input: bcf or vcf
## output: text file with tables summarizing variants and (if required) plots of the 
##         corresponding statistics
## NOTE: vcftools not running on server :( Will look in to the matter when I have the mood
#
#produce_stats(){


####################
# DEFINE VARIABLES #
####################

. /mnt/nfs/projects/ecosystem_biology/esb.shortcuts
OUTDIR=`pwd` # change?
#echo $INDIR
INPUT= 
OUTPUT='vcf'
TESTING=0
GEN=

while getopts "ht:i:o:r:" OPTION
do
    case $OPTION in
        h) usage ;;
        i) INPUT=$OPTARG ;;
        o) OUTPUT=$OPTARG ;;
        t) TESTING=$OPTARG ;;
        r) GEN=$OPTARG ;;
        ?) usage ;;
    esac
done

####################

######################
# REFERENCE SEQUENCE #
######################
## NOTE: Will change later once the metagenome is ready
if [[ "$GEN" = "cm" ]]
then
  echo "Using Candidatus microthrix genome as reference"
  REF=$MTXREF
  echo "$REF"
#elif [[ "$REF" = "ZG" ]]
elif [[ "$GEN" = "zg" ]]
then
  echo "Using Zooglea genome as reference"
  REF=$ZGLREF
  echo "$REF"
elif [[ "$GEN" = "cm_old" ]]
then
  echo "Using Zooglea OLD genome as reference"
  REF=$MTXREF_OLD
  echo "$REF"
elif [ -z "$GEN" ]
then
echo "Error: Reference genome not specified"
else
  echo "Using LAO metagenome assembly as reference"
  REF=$GEN
fi
######################

#################
# RUN FUNCTIONS #
#################
if [[ "$INPUT" =~ "bam" && "$OUTPUT" = "vcf" ]] # BAM to VCF
then
  # Run all these variant callers in parallel. Take up 3 processors when
  # running this script
run_mpileup $TESTING $REF $INPUT & 2> errormpu.log;
run_freebayes $TESTING $REF $INPUT & 2> errorfrb.log;
#run_gatk $TESTING $REF $INPUT & 2> errorgatk.log;
wait
elif [[ "$INPUT" =~ "bam" && "$OUTPUT" =~ "fasta" ]] # BAM to FASTA
then 
  echo "Input is bam"
  echo "Output is fasta"
  generate_fastq $TESTING $REF $INPUT
  FASTQ=`echo "$BAM" | sed 's/bam/fq/'` # Output file name
  generate_fasta $TESTING $FASTQ
elif [[ "$INPUT" =~ "fastq" && "$OUTPUT" =~ "fasta" ]] # FASTQ to FASTA
then 
  echo "Input is fastq"
  echo "Output is fasta"
 generate_fasta $TESTING $INPUT
else 
  # Exception
  echo "Input should be in BAM/BCF/FASTQ format"
  usage
fi
exit 0
#if [[ "$INPUT" =~ "bam" && "$OUTPUT" = "bcf" ]] # BAM to BCF
#then
#  echo "Input is bam"
#  echo "Output is bcf"
#  generate_bcf $TESTING $REF $INPUT
#elif [[ "$INPUT" =~ "bam" && "$OUTPUT" = "vcf" ]] # BAM to VCF
#then
#  echo "Input is bam"
#  echo "Output is vcf"
#  generate_bcf $TESTING $REF $INPUT
#  BCF=`echo $INPUT | rev | cut -f1 -d "/" | rev | sed 's/bam/bcf/'`
#  #BCF=`echo "$INPUT" | sed 's/bam/bcf/'` # Output file name
#  generate_vcf $TESTING $OUTDIR/$BCF
#elif [[ "$INPUT" =~ "bcf" && "$OUTPUT" =~ "vcf" ]]  # BCF to VCF
#then 
#  echo "Input is bcf"
#  echo "Output is vcf"
#  generate_vcf $TESTING $INPUT
#elif [[ "$INPUT" =~ "bam" && "$OUTPUT" =~ "fastq" ]] # BAM to FASTQ
#then 
#  echo "Input is bam"
#  echo "Output is fastq"
#  generate_fastq $TESTING $REF $INPUT

## NOTE: Commented out all the functions to generate FASTA files.

#else
#  echo "Input should be in BAM/BCF/FASTQ format"
#  usage
#fi
#################
