rm(list = ls())

suppressPackageStartupMessages({
  library(rstudioapi)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(d3heatmap)
  library(Rcpp)
  library(hciR)
  library(hciRdata)
  library(htmlwidgets)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

load("bulk_seq_tidy.rds")

#DDS object
pb_cd8_pre_dds <- deseq_from_tibble(cd8_pb_pre_count, cd8_pb_pre_sample, design = ~ case_control)

View(pb_cd8_pre_dds)

#Stailizing variance
rld_pb_cd8_pre <- rlog(pb_cd8_pre_dds)
View(rld_pb_cd8_pre)
#PCA
plot_pca(rld_pb_cd8_pre, "case_control", tooltip="id", width=700)
plot_pca(rld_pb_cd8_pre, "case_control", tooltip="id", width=700, pc = c(3,4))

#Res 
pb_cd8_pre_res <- results_all(pb_cd8_pre_dds, human100, alpha = 0.05, vs = "all", simplify = FALSE)

#Visualize how many padj < 0.05 
table(pb_cd8_pre_res[[1]]$padj < 0.05)  

#Volcano
plot_volcano(pb_cd8_pre_res[[1]], ggplot=TRUE, pvalue= -log10(sort(pb_cd8_pre_res[[1]]$padj)[30])) + ggtitle("pb_cd8_pre")

#Heatmap
x <- top_counts(pb_cd8_pre_res[[1]], rld_pb_cd8_pre, sort_fc = TRUE, padj=0.05)
plot_genes(x, "case_control", scale="row", fontsize_row= 6)

save(list = ls(), file = "bulk_seq_caseVSctrl_PrePbCD8.rds")
