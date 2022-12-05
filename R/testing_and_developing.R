library(tidyverse)
library(edgeR)
library(DESeq2)
# devtools::install_github("nicolau/multiDEG")
library(multiDEG)

load("data/phenodata.rda")
load("data/raw.exp.rda")

results <- DEG_analysis(raw.exp = raw.exp, phenodata = phenodata, treated = "INF", nontreated = "CTRL", class.column = "Class", adjust.method = "fdr")

plot(listDEGs = results, type = "up")
plot(listDEGs = results, type = "down")

