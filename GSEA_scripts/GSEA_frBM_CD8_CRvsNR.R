rm(list = ls())

#Load packages 
suppressPackageStartupMessages({
  library(DOSE)
  library(clusterProfiler)
  library(tidyverse)
  library(readr)
  library(DESeq2)
  library(fgsea)
  library(hypeR)
  library(gprofiler2)
  library(magrittr)
  library(data.table)
  library(biomaRt)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggupset)
  library(gmt)
  library(hciR)
  library(DESeq2)
  library(hciRdata)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

ts <- timestamp(stamp = Sys.Date(), frfix = "", suffix = "")
GSDirectory <- paste(PrimaryDirectory, paste(ts, "GSEA_dir", sep = "_"), sep = "/")
dir.create(GSDirectory)
load("bulk_seq_tidy.rds")

bm_cd8_fr_dds <- deseq_from_tibble(cd8_bm_fr_count, cd8_bm_fr_sample, design = ~ response)
rld_bm_cd8_fr <- rlog(bm_cd8_fr_dds)
bm_cd8_fr_res <- results_all(bm_cd8_fr_dds, human100, alpha = 0.05, vs = "all", simplify = FALSE)

#create ranks
data <- bm_cd8_fr_res[[2]]
rnk <-   data %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(gene_name, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene_name) %>%
  summarise(stat = mean(stat))
rev_rank <- rnk
ranks <- deframe(rev_rank)
ranks <- sort(ranks, decreasing = TRUE)

rnk_up <- (ranks > 0)
rnk_up <- names(rnk_up)
rnk_down <- (ranks < 0)
rnk_down <- names(rnk_down)


rnk_up_entrez <- bitr(rnk_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
rnk_down_entrez <- bitr(rnk_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#GO analysis
ego <- enrichGO(names(ranks), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

head(ego)
goplot(ego)
barplot(ego, showCategory = 30)
dotplot(ego, showCategory = 30)

#Reduce redundancy 
sego <- simplify(ego)
cnetplot(sego, foldChange = ranks)
cnetplot(sego, foldChange = ranks, circular = TRUE, colorEdge = TRUE)
upsetplot(sego)
emapplot(sego)
heatplot(sego, foldChange = ranks)

# pathways
hallmark <- gmtPathways("h.all.v7.1.symbols.gmt")

fgseaRes <- fgsea(hallmark, ranks, maxSize = 500, minSize = 15, eps = 0)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Plot enrichment
plotEnrichment(hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], ranks) + labs(title = "INTERFERON_ALPHA_RESPONSE")
plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks) + labs(title = "INTERFERON_GAMMA_RESPONSE")
plotEnrichment(hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) + labs(title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

#C7 hallmark (immune signature)
hallmark_C7 <- gmtPathways("c7.all.v7.1.symbols.gmt")

fgseaRes <- fgsea(hallmark_C7, ranks, maxSize = 500, minSize = 15, eps = 0)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#KEGG analysis
#genes_up <- deframe(rnk_up_entrez)
#genes_up <- sort(genes_up, decreasing = TRUE)
#kk <- enrichKEGG(genes_up, organism = "human")
#gk <- gseKEGG(genes_up, organism="human", nPerm=10000)

save(list = ls(), file = "GSEA_frBM_CD8_CRvsNR.rds")
