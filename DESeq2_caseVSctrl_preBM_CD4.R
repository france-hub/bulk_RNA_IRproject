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
  library(htmlwidgets)
})

load("bulk_seq_tidy.rds")

#DDS object
bm_cd4_pre_dds <- deseq_preom_tibble(cd4_bm_pre_count, cd4_bm_pre_sample, design = ~ case_control)

#Stailizing variance
rld_bm_cd4_pre <- rlog(bm_cd4_pre_dds)

#PCA
plot_pca(rld_bm_cd4_pre, "case_control", tooltip="id", width=700)
plot_pca(rld_bm_cd4_pre, "case_control", tooltip="id", width=700, pc = c(3,4))

#Res 
bm_cd4_pre_res <- results_all(bm_cd4_pre_dds, human100, alpha = 0.05, vs = "all", simplify = FALSE)

#Visualize how many padj < 0.05 
table(bm_cd4_pre_res[[1]]$padj < 0.05)  

#Volcano
plot_volcano(bm_cd4_pre_res[[1]], ggplot=TRUE, pvalue= -log10(sort(bm_cd4_pre_res[[1]]$padj)[30])) + ggtitle("bm_cd4_pre")

#Heatmap
x <- top_counts(bm_cd4_pre_res[[1]], rld_bm_cd4_pre, sort_fc = TRUE, padj=0.05)
plot_genes(x, "case_control", scale="row", fontsize_row= 6)

save(list = ls(), file = "bulk_seq_caseVSctrl_prebmcd4.rds")
