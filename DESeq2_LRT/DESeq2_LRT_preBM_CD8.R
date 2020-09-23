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
bm_cd8_pre_dds <- deseq_from_tibble(cd8_bm_pre_count, cd8_bm_pre_sample, design = ~ response)

#Create DESeq object
bm_cd8_pre_LRT <- DESeq(bm_cd8_pre_dds, test="LRT", reduced = ~ 1)

#Stailizing variance
rld_Prebm_cd8_LRT <- rlog(bm_cd8_pre_LRT)

#Results
res_Prebm_cd8_LRT <- results(bm_cd8_pre_LRT)

##Filter only genes which are significant
padj.cutoff <-  0.05
res_Prebm_cd8_LRT <- res_Prebm_cd8_LRT %>%
  data.frame() %>%
  filter(padj < padj.cutoff)

##Add columns with SYMBOL (adapt tibble for hciR package)
biomart <- read_biomart("human")
Columns <- c("gene_name", "biotype", "chromosome",  "description")
if("human_homolog" %in% names(biomart)) Columns <- c(Columns, "human_homolog")
res_Prebm_cd8_LRT <- annotate_results(res_Prebm_cd8_LRT, biomart, Columns)

x <- top_counts(res_Prebm_cd8_LRT, rld_Prebm_cd8_LRT, sort_fc = FALSE, filter = FALSE, padj=0.05)
plot_genes(x, "response", scale="row", fontsize_row= 6)

res_Prebm_cd8_LRT %>% 
  arrange(padj) %>% 
  DT::datatable()

save(list = ls(), file = "bulk_seq_LRT_Prebmcd8.rds")
load("bulk_seq_LRT_Prebmcd8.rds")
