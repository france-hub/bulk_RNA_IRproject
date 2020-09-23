rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(d3heatmap)
  library(Rcpp)
  library(hciR)
  library(hciRdata)
  library(biomaRt)
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
pb_cd8_pre_dds <- deseq_from_tibble(cd8_pb_pre_count, cd8_pb_pre_sample, design = ~ response)

#Create DESeq object
pb_cd8_pre_LRT <- DESeq(pb_cd8_pre_dds, test="LRT", reduced = ~ 1)

#Stailizing variance
rld_Prepb_cd8_LRT <- rlog(pb_cd8_pre_LRT)

#Results
res_Prepb_cd8_LRT <- results(pb_cd8_pre_LRT)

##Filter only genes which are significant
padj.cutoff <-  0.05
res_Prepb_cd8_LRT <- res_Prepb_cd8_LRT %>%
  data.frame() %>%
  filter(padj < padj.cutoff)

##Add columns with SYMBOL (adapt tibble for hciR package)
biomart <- read_biomart("human")
Columns <- c("gene_name", "biotype", "chromosome",  "description")
if("human_homolog" %in% names(biomart)) Columns <- c(Columns, "human_homolog")
res_Prepb_cd8_LRT <- annotate_results(res_Prepb_cd8_LRT, biomart, Columns)

x <- top_counts(res_Prepb_cd8_LRT, rld_Prepb_cd8_LRT, sort_fc = FALSE, filter = FALSE, padj=0.05)
plot_genes(x, "response", scale="row", fontsize_row= 6)

#Volcano
res_volc <- results_all(pb_cd8_pre_dds, human100, alpha = 0.05, vs = "all", simplify = FALSE)
plot_volcano(res_volc[[3]], ggplot=TRUE, pvalue= -log10( sort(res_volc[[3]]$padj)[50])) + ggtitle("pb_cd8_pre")


res_table <- res_Prepb_cd8_LRT %>% 
  arrange(padj) %>% 
  DT::datatable()

save(list = ls(), file = "bulk_seq_LRT_PrePbCD8.rds")
load("bulk_seq_LRT_PrePbCD8.rds")
