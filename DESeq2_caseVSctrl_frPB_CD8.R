rm(list = ls())

supfrssPackageStartupMessages({
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
pb_cd8_fr_dds <- deseq_from_tibble(cd8_pb_fr_count, cd8_pb_fr_sample, design = ~ case_control)

#Stailizing variance
rld_pb_cd8_fr <- rlog(pb_cd8_fr_dds)

#PCA
plot_pca(rld_pb_cd8_fr, "case_control", tooltip="id", width=700)
plot_pca(rld_pb_cd8_fr, "case_control", tooltip="id", width=700, pc = c(3,4))

#Res 
pb_cd8_fr_res <- results_all(pb_cd8_fr_dds, human100, alpha = 0.05, vs = "all", simplify = FALSE)

#Visualize how many padj < 0.05 
table(pb_cd8_fr_res[[1]]$padj < 0.05)  

#Volcano
plot_volcano(pb_cd8_fr_res[[1]], ggplot=TRUE, pvalue= -log10(sort(pb_cd8_fr_res[[1]]$padj)[30])) + ggtitle("pb_cd8_fr")

#Heatmap
x <- top_counts(pb_cd8_fr_res[[1]], rld_pb_cd8_fr, sort_fc = TRUE, padj=0.05)
plot_genes(x, "case_control", scale="row", fontsize_row= 6)

save(list = ls(), file = "bulk_seq_caseVSctrl_frPbCD8.rds")
