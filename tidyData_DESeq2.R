rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################
if(!require('DESeq2')) {
  install.packages('DESeq2')
}
if(!require('readxl')){
  install.packages('readxl')
}
if(!require('devtools')){
  install.packages('devtools')
}
if(!require('ririzarr/rafalib')){
  devtools::install_github('ririzarr/rafalib')
}
if(!require('vsn')){
  install.packages('vsn')
}
if(!require('ggplot2')){
  install.packages('ggplot2')
}
if(!require('pheatmap')){
  install.packages('pheatmap')
}
if(!require('org.Hs.eg.db')){
  BiocManager::install("org.Hs.eg.db")
}
if(!require('org.Hs.eg.db')){
  BiocManager::install("org.Hs.eg.db")
}
if(!require('dplyr')){
  BiocManager::install("dplyr")
}
if(!require('tidyverse')){
  BiocManager::install("tidyverse")
}
if(!require('EnhancedVolcano')){
  BiocManager::install("EnhancedVolcano")
}
if(!require("rstudio/d3heatmap")){
  devtools::install_github("rstudio/d3heatmap")
}
if(!require("Rcpp")){
  install.packages("Rcpp")
}
if(!require("HuntsmanCancerInstitute/hciR")){
  devtools::install_github("HuntsmanCancerInstitute/hciR")
}
if(!require("HuntsmanCancerInstitute/hciRdata")){
  devtools::install_github("HuntsmanCancerInstitute/hciRdata")
}
###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(devtools)
  library(rafalib)
  library(vsn)
  library(ggplot2)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(hgu95av2.db)
  library(dplyr)
  library(EnhancedVolcano)
  library(tidyverse)
  library(d3heatmap)
  library(Rcpp)
  library(hciR)
  library(hciRdata)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define workingDirectory
wdName <- "work_dir"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

# Counts_metadata dir
cm_name <- "counts_metadata"
counts_metadata <- paste(PrimaryDirectory, cm_name, sep = "/")
dir.create(counts_metadata)

#Read sample metadata and change complete in CR and no in NR
samples <- read_tsv(paste(counts_metadata, "samples.txt", sep="/"))
samples$response <- gsub("complete", "CR", samples$response)
samples$response <- gsub("no", "NR", samples$response)

samples$response <- as.factor(samples$response)
levels(samples$response)

#Read counts
count <- read_tsv(paste(counts_metadata, "counts.txt", sep="/"))

#plot_filter(counts) 

#Change column names
colnames(count) <-  gsub("+", "", colnames(count), fixed = TRUE)
colnames(count) <-  gsub("fr1", "fr", colnames(count))
colnames(count) <-  gsub("fr2", "fr", colnames(count))
colnames(count)[duplicated(colnames(count))] <- gsub("fr", "fr2", colnames(count)[duplicated(colnames(count))])

# change sample$id names
samples$id <-  gsub("fr1", "fr", samples$id)
samples$id <-  gsub("fr2", "fr", samples$id)
samples$id[duplicated(samples$id)] <- gsub("fr", "fr2", samples$id[duplicated(samples$id)])
samples$id %in% colnames(count) #check correspondence

#concordance sample metadata and count matrix
# check for concordance between samples and counts data frames
x <-  sort(samples$id)
y <-  sort(colnames(count)[2:ncol(count)])

xy <-  as.data.frame(cbind(x,y))

table(x == y)

xy[which(xy$x != xy$y) ,]

# expected dimensions
table(samples$cell_type, samples$organ)

##Subset different counts per timepoint and cell type

#Pre
pre_pb_cd4 <- samples$id[grepl("pbprecd4|hdpb[0-9]cd4|hdpb[a-z]{2}cd4", samples$id)]
cd4_pb_pre_count <-  count[, c("geneid", pre_pb_cd4)]

pre_pb_cd8 <- samples$id[grepl("pbprecd8|hdpb[0-9]cd8|hdpb[a-z]{2}cd8", samples$id)]
cd8_pb_pre_count <-  count[, c("geneid", pre_pb_cd8)]

pre_bm_cd4 <- samples$id[grepl("bmprecd4|hdbm[0-9]{4}cd4|hd[a-z]{2}bmcd4", samples$id)]
cd4_bm_pre_count <-  count[, c("geneid", pre_bm_cd4)]

pre_bm_cd8 <- samples$id[grepl("bmprecd8|hdbm[0-9]{4}cd8|hd[a-z]{2}bmcd8", samples$id)]
cd8_bm_pre_count <-  count[, c("geneid", pre_bm_cd8)]

#FR
fr_pb_cd4 <- samples$id[grepl("pbfrcd4|hdpb[0-9]cd4|hdpb[a-z]{2}cd4", samples$id)]
cd4_pb_fr_count <-  count[, c("geneid", fr_pb_cd4)]

fr_pb_cd8 <- samples$id[grepl("pbfrcd8|hdpb[0-9]cd8|hdpb[a-z]{2}cd8", samples$id)]
cd8_pb_fr_count <-  count[, c("geneid", fr_pb_cd8)]

fr_bm_cd4 <- samples$id[grepl("bmfrcd4|hdbm[0-9]{4}cd4|hd[a-z]{2}bmcd4", samples$id)]
cd4_bm_fr_count <-  count[, c("geneid", fr_bm_cd4)]

fr_bm_cd8 <- samples$id[grepl("bmfrcd8|hdbm[0-9]{4}cd8|hd[a-z]{2}bmcd8", samples$id)]
cd8_bm_fr_count <-  count[, c("geneid", fr_bm_cd8)]

##Subsets samples
#Pre
cd4_pb_pre_sample <- samples[samples$id %in% pre_pb_cd4,]
cd8_pb_pre_sample <- samples[samples$id %in% pre_pb_cd8,]
cd4_bm_pre_sample <- samples[samples$id %in% pre_bm_cd4,]
cd8_bm_pre_sample <- samples[samples$id %in% pre_bm_cd8,]

#FR
cd4_pb_fr_sample <- samples[samples$id %in% fr_pb_cd4,]
cd8_pb_fr_sample <- samples[samples$id %in% fr_pb_cd8,]
cd4_bm_fr_sample <- samples[samples$id %in% fr_bm_cd4,]
cd8_bm_fr_sample <- samples[samples$id %in% fr_bm_cd8,]

save(list = ls(), file = "bulk_seq_tidy.rds")

