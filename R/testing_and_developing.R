library(tidyverse)
library(edgeR)
library(DESeq2)
library(multiDEG)

load("data/phenodata.rda")
load("data/raw.exp.rda")

results <- DEG_analysis(raw.exp = raw.exp, phenodata = phenodata, treated = "INF", nontreated = "CTRL", class.column = "Class", adjust.method = "fdr")


genes_wilcox_up   <- results$`Wilcoxon rank-sum test` %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(log2foldChange > 1  & pValues < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)
genes_wilcox_down <- results$`Wilcoxon rank-sum test` %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(log2foldChange < -1 & pValues < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)


genes_edgeR_up   <- results$edgeR %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(logFC > 1  & PValue < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)
genes_edgeR_down <- results$edgeR %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(logFC < -1 & PValue < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)


genes_DESeq2_up   <- results$DESeq2 %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(log2FoldChange > 1  & pvalue < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)
genes_DESeq2_down <- results$DESeq2 %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::filter(log2FoldChange < -1 & pvalue < 0.05) %>%
  dplyr::select(Symbol) %>%
  unlist(use.names = F)

source("https://raw.githubusercontent.com/nicolau/code-R/master/funcoes_para_diagrama_venn.R")

plot.triple.venn(a1 = genes_wilcox_down, a2 = genes_edgeR_down, a3 = genes_DESeq2_down,
                 labels = c("Wilcox", "edgeR", "DESeq2"))

plot.triple.venn(a1 = genes_wilcox_up, a2 = genes_edgeR_up, a3 = genes_DESeq2_up,
                 labels = c("Wilcox", "edgeR", "DESeq2"))
