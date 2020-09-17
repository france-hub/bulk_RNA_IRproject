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

load("bulk_seq_tidy.rds")

#DDS object
pb_cd4_pre_dds <- deseq_from_tibble(cd4_pb_pre_count, cd4_pb_pre_sample, design = ~ response)

#Create DESeq object
pb_cd4_pre_LRT <- DESeq(pb_cd4_pre_dds, test="LRT", reduced = ~ 1)

#Stailizing variance
rld_Prepb_cd4_LRT <- rlog(pb_cd4_pre_LRT)

#Results
res_Prepb_cd4_LRT <- results(pb_cd4_pre_LRT)

##Filter only genes which are significant
padj.cutoff <-  0.05
res_LRT <- res_Prepb_cd4_LRT %>%
  data.frame() %>%
  filter(padj < padj.cutoff)

##Add columns with SYMBOL (adapt tibble for hciR package)
biomart <- read_biomart("human")
Columns <- c("gene_name", "biotype", "chromosome",  "description")
if("human_homolog" %in% names(biomart)) Columns <- c(Columns, "human_homolog")
res_LRT <- annotate_results(res_LRT, biomart, Columns)

x <- top_counts(res_LRT, rld_Prepb_cd4_LRT, sort_fc = FALSE, filter = FALSE, padj=0.05)
plot_genes(x, "response", scale="row", fontsize_row= 6)

res_table <- res_LRT %>% 
  arrange(padj) %>% 
  DT::datatable()

save(list = ls(), file = "bulk_seq_LRT_PrePbcd4.rds")
load("bulk_seq_LRT_PrePbcd4.rds")
