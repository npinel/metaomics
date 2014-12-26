#!/bin/R

##########################################################################################
# Author(s):                                                                             #
#    1. Shaman Narayanasamy (shaman.narayanasamy@uni.lu)                                 #
# Affiliation(s):                                                                        #
#    1. Luxembourg Center for Systems Biomedicine                                        #
# Project(s):                                                                            #
#    Time-resolved multi-omic analysis of lipid-accumulating communities                 #
# Date: 03 November 2013                                                                 #
# Script: 04b.analysis.idba.sh                                                           #
# Purpose:                                                                               #
#    R script for reading in cufflinks data sets                                         #
# NOTES:                                                                                 #
#    TO DO: 1. Make script more dynamic. Provide arguements for the variables.           #
#           2. Rearrange the time points accordingly and match RNA and DNA measurements  #
#           3. Use R packages for analysis? Or implement custom analysis (discuss with   #
#              Nic)                                                                      #
#           4. Make a wrapper shell script                                               # 
##########################################################################################

##########################################################################################
# Function: Read in all the cufflinks fpkm values based on a list of files and merge 
# the results into a single time series table.
# NOTE: Same function can be used for MG or MT data sets provided the mappings and 
#       corresponding annotation file is used
##########################################################################################

ts.read <- function(file.list){

  M <- read.delim(file.list[1])[,c("locus","FPKM")] 
    for (i in 2:length(file.list)){
	M1 <- read.delim(file.list[i])[,c("locus","FPKM")]
	M <- merge(M,M1,by="locus")
    }

  colnames(M)[2:ncol(M)] <- file.list[1:length(file.list)]
return(M)
}

##########################################################################################
# Function: Integrate MT and MG data sets. Performs RNA/DNA normalization for whatever
#	    sets of samples provided resulting in a single normalized matrix.
# NOTE: Draft version.
# TODO: 1. Data transformation (log)
#	2. Handling 0 expression (MT) values (division of 0)
#       3. Handling 0 MG values (division by 0)
##########################################################################################

MTMG.normalize <- function(MT.mat,MG.mat){

# Simplify sample names MG
colnames(MG.mat)[2:ncol(MG.mat)] <- gsub(pattern="_MG.genes.fpkm_tracking",replacement="",x=colnames(MG.mat)[2:ncol(MG.mat)])

# Simplify sample names MG
colnames(MT.mat)[2:ncol(MT.mat)] <- gsub(pattern="_MT.genes.fpkm_tracking",replacement="",x=colnames(MT.mat)[2:ncol(MT.mat)])

# Find intersection between columns of MT and MG data to make sure the samples are paired
MGMT.paired <- intersect(colnames(MT.mat)[2:ncol(MT.mat)],colnames(MG.mat)[2:ncol(MG.mat)])

# Normalize the matrix by RNA/DNA measure
MTMG.norm.mat<- cbind(MT.mat[,1],MT.mat[,MGMT.paired]/MG.mat[,MGMT.paired])

colnames(MTMG.norm.mat)[1] <- "locus"

return(MTMG.norm.mat)
}

##########################################################################################
# Function: Integrate MT and MG data sets. Performs RNA/DNA normalization for whatever
#	    sets of samples provided resulting in a single normalized matrix.
# NOTE: Draft version.
# TODO: 1. Data transformation (log)
#	2. Handling 0 expression (MT) values (division of 0)
#       3. Handling 0 MG values (division by 0)
##########################################################################################
require(stringr)

gff.unlist=function(gff.file){
x <- read.delim(gff.file,header=F)
y <- cbind(str_split_fixed(x[,dim(x)[2]],";",3))
y[,1] <- gsub("ID=","",y[,1])
y[,2] <- gsub("Name=","",y[,2])
y[,11] <- gsub("Ontology_term=","",y[,11])
y[which(annot[,2]==""),2]=NA
y[which(annot[,3]==""),3]=NA
y[which(annot[,10]==""),10]=NA
y[which(annot[,11]==""),11]=NA
x <- cbind(x[,-dim(x)[2]],y)
x<-cbind(x,str_split_fixed(x[,1],":",2)[,1])
colnames(x) <- c("ID","source","feature","start","stop","score","strand","frame","feature_ID","function(EC)","Ontology_term","contig")
x[,1] <- paste(x[,1],":",x[,4]-1,"-",x[,5],sep="")
return(x)
}

##########################################################################################
# Function: Handling exceptions in the annotation of transcripts
#
# NOTE: This came about when I realized that cufflinks seems to call transcripts that 
#       were not in the annotation file (gff). In most these cases, cufflinks concatenates
#       two transcripts that are close by (loci), while those loci, do not exist within
#       the gff file. Cufflinks also calls transcripts that do not even exist at all in
#       the annotation file. For now, we can just omit those files.
# TODO: 1. Clean up
#       2. Create a table with all the "bad ids"
#	3. Convert to function
#       4. Figure out what to do with the missing values
##########################################################################################

# create a test matrix
test <- matrix(nrow=nrow(bad.ids.3),ncol=3)
# start loop
for (i in 1:nrow(bad.ids.3)){
# extract the contig containing a given transcript
  sub <- annot.2[which(annot.2$contig==as.character(bad.ids.3$contig[i])),]
  print(length(which(sub$start==bad.ids.3$start[i])))
  start.loc <- length(which(sub$start==bad.ids.3$start[i]))
  print(length(which(sub$stop==bad.ids.3$stop[i])))
  stop.loc <- length(which(sub$stop==bad.ids.3$stop[i]))
  if (start.loc > 0 & stop.loc > 0){
    id.1 <- as.character(sub[which(sub$start==bad.ids.3$start[i]),"feature_ID"])
    id.2 <- as.character(sub[which(sub$stop==bad.ids.3$stop[i]),"feature_ID"])
    test[i,] <- c(as.character(bad.ids.3$cufflink_id[i]),id.1,id.2)
  }
  else{
    test[i,] <- c(as.character(bad.ids.3$cufflink_id[i]),NA,NA)
  }
print(test[i,])
}
bad.ids.4 <- cbind(bad.ids.3,test[,2:3])

colnames(bad.ids.4)[5:6] <- c("ID_1","ID_2")


bad.ids.5<-cbind(as.character(bad.ids.4$cufflink_id),NA,NA,as.character(bad.ids.4$start),as.character(bad.ids.4$stop),NA,NA,NA,paste(bad.ids.4$ID_1,bad.ids.4$ID_2,sep=";"),NA,NA,as.character(bad.ids.4$contig))
colnames(bad.ids.5) <- colnames(annot.2)
annot.3<-rbind(annot.2,bad.ids.5)

##########################################################################################
# Function: Main
# NOTE: At the moment, this is all the main script does. Possible to provide more options
#       depending on the input.
# TODO: Provide more options for expression matrices.
#	1. Reading in a paired sample data set
#	2. Reading in a time series data set
##########################################################################################

args <- commandArgs(trailingOnly=TRUE)
print("Reading in list of fpkm_tracking files: args[1]")_
system("date")
file.list <- as.character(read.table(args[1],header=FALSE)[,1])

MT.mat <- ts.read(MT.file.list)
MG.mat <- ts.read(MG.file.list)
E.mat <- MTMG.normalize(MT.mat,MG.mat)

print("Writing out expression matrix")
out.file <- args[2]
write.table(E.mat,out.file,quote=F,sep="\t",row.names=F)
system("date")

MG.cols<-str_split_fixed(colnames(MG.mat),"/",11)[-1,11]
MG.cols<-gsub("_MG.genes.fpkm_tracking","",MG.cols)
MT.cols<-str_split_fixed(colnames(MT.mat),"/",11)[-1,11]
MT.cols<-gsub("_MT.genes.fpkm_tracking","",MT.cols)




