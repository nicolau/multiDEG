#' Calculates DE genes table for Wilcoxon rank-sum test, DESeq2 and edgeR
#'
#' @param raw.exp data.frame of raw counts
#' @param phenodata data.frame of phenotypic data for all samples
#' @param treated A string containing the name of class to be considered as treatment
#' @param nontreated A string containing the name of class to be considered as nontreatment
#' @param class.column A string contatining the name of column where treated and nontreated will be extracted
#' @param adjust.method Merhod of multiple comparison will be used for Wilcoxon rank-sum test pvalues.
#' @param covariables A vector list of variables should be used to covariates adjusment.
#' @param paired.samples A boolean variable to indicate paired-sample comparison
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
DEG_analysis <- function(raw.exp, phenodata, treated, nontreated, class.column = "Class", adjust.method = "fdr", covariables = NULL, paired.samples.column = NULL) {

  # class.column          <- "Class"
  # adjust.method         <- "fdr"
  # treated               <- "INF"
  # nontreated            <- "CTRL"
  # covariables           <- c("Gender")
  # covariables           <- NULL
  # paired.samples.column <- NULL
  # data(phenodata)
  # data(raw.exp)
  #
  # library(tidyverse)
  # library(DESeq2)
  # library(edgeR)


  covariablesStop   <- FALSE
  covariablesVector <- NULL
  for(var in covariables) {
    if(!var %in% colnames(phenodata)) {
      covariablesStop <- TRUE
      covariablesVector <- c(covariablesVector, var)
    }
  }
  if(covariablesStop) {stop("Covariables: \"", paste(covariablesVector, collapse = ", "), "\" has not found in the phenodata file!")}

  result <- list()

  if(class.column != "Class") {
    if(sum(colnames(phenodata) == "Class") == 1) {
      colnames(phenodata)[colnames(phenodata) == "Class"] <- "ClassOld"
    }
    colnames(phenodata)[colnames(phenodata) == class.column] <- "Class"
  }


  if(length(covariables) > 1) {
    phenodata <- phenodata magrittr::`%>%` tidyr::unite(Concatenate, c(covariables), remove = FALSE)
    model <- paste0("~0+", paste(c(paired.samples.column, "Concatenate", class.column), collapse = "+"))
  } else if(length(covariables) == 1) {
    model <- paste0("~0+", paste(c(paired.samples.column, covariables, class.column), collapse = "+"))
  } else {
    model <- paste0("~0+", paste(c(paired.samples.column, class.column), collapse = "+"))
  }


  phenodata  <- phenodata magrittr::`%>%` dplyr::filter(Class %in% c(treated, nontreated))
  raw.exp    <- raw.exp magrittr::`%>%` dplyr::select(phenodata$Sample)

  conditions <- factor(t(phenodata$Class)) magrittr::`%>%` relevel(conditions, ref = nontreated)

  # Use this variable as covariates
  # gender     <- factor(t(phenodata$Gender))

  if(!identical(colnames(raw.exp), phenodata$Sample)) {
    stop("Error: colnames(raw.exp) should be the same of phenodata$Sample")
  }

  ############################################################ Wilcoxon rank-sum test ############################################################
  message("Calculating DE genes using Wilcoxon rank-sum test...")
  # edgeR TMM normalize
  y <- edgeR::DGEList(counts = raw.exp, group = conditions)
  ## Remove rows conssitently have zero or very low counts
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, keep.lib.sizes = FALSE]
  ## Perform TMM normalization and convert to CPM (Counts Per Million)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  norm.exp <- edgeR::cpm(y) magrittr::`%>%` as.data.frame()
  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(norm.exp), function(i){
    data <- cbind.data.frame(gene = as.numeric(t(norm.exp[i,])), conditions)
    p <- wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = adjust.method)
  # Calculate the fold-change for each gene
  dataCon1 <- norm.exp magrittr::`%>%` dplyr::select(c(which(conditions==nontreated)))
  dataCon2 <- norm.exp magrittr::`%>%` dplyr::select(c(which(conditions==treated)))
  foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  # Output results based on the FDR threshold 0.05
  outRst <- data.frame(row.names = rownames(norm.exp), log2FoldChange = foldChanges, pvalue = pvalues, padj = fdr) magrittr::`%>%` na.omit(outRst)

  result[['Wilcoxon rank-sum test']] <- outRst
  message("Done!")
  message("")
  ############################################################ Wilcoxon rank-sum test ############################################################


  #################################################################### DESeq2 ####################################################################
  message("Calculating DE genes using DESeq2...")
  contrasts <- c("Class", treated, nontreated)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw.exp, colData = phenodata, design = ~ Class)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw.exp, colData = phenodata, design = as.formula(model))

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrasts)

  tTags <- res magrittr::`%>%` as.data.frame()
  result[['DESeq2']] <- tTags
  message("Done!")
  message("")
  #################################################################### DESeq2 ####################################################################


  ##################################################################### edgeR ####################################################################
  message("Calculating DE genes using edgeR...")
  group  <- conditions
  y      <- edgeR::DGEList(raw.exp)
  y      <- edgeR::calcNormFactors(y)
  term   <- ~0+group
  # term   <- as.formula(model)
  design <- model.matrix(term)

  y      <- edgeR::estimateDisp(y, design)
  fit    <- edgeR::glmFit(y, design)
  lrt    <- edgeR::glmLRT(fit, coef = 2, contrast = c(-1, 1))
  tTags  <- edgeR::topTags(lrt, n = NULL)[["table"]]
  result[['edgeR']] <- tTags
  message("Done!")
  message("")
  ##################################################################### edgeR ####################################################################


  ################################################################# Overlap DEGs #################################################################
  overlap_result <- overlap_DEGs(listDEGs = result)
  result[["overlap"]] <- overlap_result
  ################################################################# Overlap DEGs #################################################################

  result[['edgeR']] <- result[['edgeR']] magrittr::`%>%` dplyr::rename(log2FoldChange = logFC, pvalue = PValue, padj = FDR)

  return(result)
}

#' Calculates overlapping genes between Wilcox rank-sum test, DESeq2 and edgeR DE methods
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
overlap_DEGs <- function(listDEGs, p_cutoff = 0.05, log2fc_cutoff = 1, padjusted = F) {

  genes_wilcox_up   <- listDEGs$`Wilcoxon rank-sum test` magrittr::`%>%` tibble::rownames_to_column("Symbol")
  genes_wilcox_down <- listDEGs$`Wilcoxon rank-sum test` magrittr::`%>%` tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_wilcox_up <- genes_wilcox_up magrittr::`%>%`
      dplyr::filter(log2foldChange > log2fc_cutoff  & FDR < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)

    genes_wilcox_down <- genes_wilcox_down magrittr::`%>%`
      dplyr::filter(log2foldChange < -log2fc_cutoff & FDR < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  } else {
    genes_wilcox_up <- genes_wilcox_up magrittr::`%>%`
      dplyr::filter(log2foldChange > log2fc_cutoff  & pValues < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)

    genes_wilcox_down <- genes_wilcox_down magrittr::`%>%`
      dplyr::filter(log2foldChange < -log2fc_cutoff & pValues < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  }


  genes_edgeR_up   <- listDEGs$edgeR magrittr::`%>%` tibble::rownames_to_column("Symbol")
  genes_edgeR_down <- listDEGs$edgeR magrittr::`%>%` tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_edgeR_up<- genes_edgeR_up magrittr::`%>%`
      dplyr::filter(logFC > log2fc_cutoff  & FDR < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
    genes_edgeR_down <- genes_edgeR_down magrittr::`%>%`
      dplyr::filter(logFC < -log2fc_cutoff & FDR < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  } else {
    genes_edgeR_up<- genes_edgeR_up magrittr::`%>%`
      dplyr::filter(logFC > log2fc_cutoff  & PValue < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
    genes_edgeR_down <- genes_edgeR_down magrittr::`%>%`
      dplyr::filter(logFC < -log2fc_cutoff & PValue < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  }


  genes_DESeq2_up   <- listDEGs$DESeq2 magrittr::`%>%` tibble::rownames_to_column("Symbol")
  genes_DESeq2_down <- listDEGs$DESeq2 magrittr::`%>%` tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_DESeq2_up <- genes_DESeq2_up magrittr::`%>%`
      dplyr::filter(log2FoldChange > log2fc_cutoff  & padj < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
    genes_DESeq2_down <- genes_DESeq2_down magrittr::`%>%`
      dplyr::filter(log2FoldChange < -log2fc_cutoff & padj < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  } else {
    genes_DESeq2_up <- genes_DESeq2_up magrittr::`%>%`
      dplyr::filter(log2FoldChange > log2fc_cutoff  & pvalue < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
    genes_DESeq2_down <- genes_DESeq2_down magrittr::`%>%`
      dplyr::filter(log2FoldChange < -log2fc_cutoff & pvalue < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  }


  plot_down <- plot.triple.venn(a1 = genes_wilcox_down, a2 = genes_edgeR_down, a3 = genes_DESeq2_down, labels = c("Wilcox", "edgeR", "DESeq2"))
  plot_up   <- plot.triple.venn(a1 = genes_wilcox_up,   a2 = genes_edgeR_up,   a3 = genes_DESeq2_up,   labels = c("Wilcox", "edgeR", "DESeq2"))

  if(!is.null(dev.list())) { dev.off() }

  return(list("plot_up" = plot_up, "plot_down" = plot_down))

}


#' Plot triple venn diagram
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
plot.triple.venn <- function(a1 = c( "a", "b", "c", "d", "e", "t", "g", "h" ),
                             a2 = c( "x", "b", "c", "d", "e", "f", "g", "z" ),
                             a3 = c( "y", "b", "c", "d", "x", "f", "g", "h" ),
                             labels = c("Group1", "Group2", "Group3")) {
  result <- list()

  exclusive <- NULL
  data      <- NULL

  data$values <- unique(c(a1, a2, a3))

  data$a1   <- a1
  data$a2   <- a2
  data$a3   <- a3

  data$n12  <- intersect(a1, a2)
  data$n13  <- intersect(a1, a3)

  data$n23  <- intersect(a2, a3)

  data$n123 <- intersect(data$n12, a3)

  # Reference four-set diagram
  venn.plot <- VennDiagram::draw.triple.venn(
    area1    = length(data$a1),
    area2    = length(data$a2),
    area3    = length(data$a3),
    n12      = length(data$n12),
    n13      = length(data$n13),
    n23      = length(data$n23),
    n123     = length(data$n123),
    category = labels,
    fill     = c("green", "red", "blue"),
    lty      = "dashed",
    cex      = 2,
    cat.cex  = 2,
    cat.col  = c("green", "red", "blue")
  )


  exclusive$a1 <- c(setdiff(setdiff(data$a1, data$a2), data$a3))
  exclusive$a2 <- c(setdiff(setdiff(data$a2, data$a1), data$a3))
  exclusive$a3 <- c(setdiff(setdiff(data$a3, data$a2), data$a1))

  data$n12 <- setdiff(data$n12, data$n123)
  data$n13 <- setdiff(data$n13, data$n123)
  data$n23 <- setdiff(data$n23, data$n123)


  result[["plot"]]            <- venn.plot
  result[["shared_all"]]      <- data$n123
  result[["shared_G1_vs_G2"]] <- data$n12
  result[["shared_G1_vs_G3"]] <- data$n13
  result[["shared_G2_vs_G3"]] <- data$n23
  result[["exclusive_G1"]]    <- exclusive$a1
  result[["exclusive_G2"]]    <- exclusive$a2
  result[["exclusive_G3"]]    <- exclusive$a3

  return(result)
}

#' Plot function
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
plot <- function(listDEGs, type = c("up", "down")) {
  if(!is.null(dev.list())) { dev.off() }
  if(type == "up") {
    grid.draw(listDEGs$overlap$plot_up$plot)
  } else {
    grid.draw(listDEGs$overlap$plot_down$plot)
  }
}

#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
get_DEG_table <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR")) {
  table <- NULL
  if(method == "Wilcox") {
    table <- results$`Wilcoxon rank-sum test`
  } else {
    table <- results[[method]]
  }
  return(table)
}

#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
get_DEG_symbols <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), direction = c("up", "down"),  p_cutoff = 0.05, log2fc_cutoff = 1, padjusted = F) {

  symbols    <- NULL
  table_DEGs <- NULL

  if(method == "Wilcox") {
    table_DEGs <- listDEGs$`Wilcoxon rank-sum test`
  } else if(method == "DESeq2") {
    table_DEGs <- listDEGs$DESeq2
  } else {
    table_DEGs <- listDEGs$edgeR
  }

  if(direction == "up") {
    symbols <- table_DEGs magrittr::`%>%`
      dplyr::filter(log2FoldChange > log2fc_cutoff) magrittr::`%>%`
      tibble::rownames_to_column("Symbol")
  } else {
    symbols <- table_DEGs magrittr::`%>%`
      dplyr::filter(log2FoldChange < -log2fc_cutoff) magrittr::`%>%`
      tibble::rownames_to_column("Symbol")
  }

  if(padjusted) {
    symbols <- symbols magrittr::`%>%`
      dplyr::filter(padj < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  } else {
    symbols <- symbols magrittr::`%>%`
      dplyr::filter(pvalue < p_cutoff) magrittr::`%>%`
      dplyr::select(Symbol) magrittr::`%>%`
      unlist(use.names = F)
  }

  return(table)
}

#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
pvalue_distribution_plot <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), padjusted = F) {
  if(method == "wilcox") {
    table_DEGs <- listDEGs$`Wilcoxon rank-sum test`
  } else if(method == "DESeq2") {
    table_DEGs <- listDEGs$DESeq2
  } else {
    table_DEGs <- listDEGs$edgeR
  }

  if(padjusted) {
    p <- ggplot2::ggplot(table_DEGs, aes(padj)) +
      ggplot2::geom_histogram() +
      ggplot2::theme_classic()
  } else {
    p <- ggplot2::ggplot(table_DEGs, aes(pvalue)) +
      ggplot2::geom_histogram() +
      ggplot2::theme_classic()
  }

  return(p)

}


#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
log2fc_correlation_plot <- function(listDEGs, method_one = c("Wilcox", "DESeq2", "edgeR"), method_two = c("Wilcox", "DESeq2", "edgeR")) {

  if(method_one != method_two) {
    if(method_one == "wilcox") {
      table_one <- listDEGs$`Wilcoxon rank-sum test` magrittr::`%>%` tibble::rownames_to_column("Symbol")
      xlabel <- "Wilcox"
    } else if(method_one == "DESeq2") {
      table_one <- listDEGs$DESeq2 magrittr::`%>%` tibble::rownames_to_column("Symbol")
      xlabel <- "DESeq2"
    } else {
      table_one <- listDEGs$edgeR magrittr::`%>%` tibble::rownames_to_column("Symbol")
      xlabel <- "edgeR"
    }

    if(method_two == "wilcox") {
      table_two <- listDEGs$`Wilcoxon rank-sum test` magrittr::`%>%` tibble::rownames_to_column("Symbol")
      ylabel <- "Wilcox"
    } else if(method_two == "DESeq2") {
      table_two <- listDEGs$DESeq2 magrittr::`%>%` tibble::rownames_to_column("Symbol")
      ylabel <- "DESeq2"
    } else {
      table_two <- listDEGs$edgeR magrittr::`%>%` tibble::rownames_to_column("Symbol")
      ylabel <- "edgeR"
    }

    table_merge <- dplyr::inner_join(x = table_one, y = table_two, by = "Symbol")

  }

  p <- ggplot2::ggplot(table_merge, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    xlab(xlabel) +
    ylab(ylabel)

  return(p)

}


#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
save_ranked_table <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), padjusted = FALSE) {

  table_DEGs <- NULL
  if(method == "wilcox") {
    table_DEGs <- listDEGs$`Wilcoxon rank-sum test` magrittr::`%>%` tibble::rownames_to_column("Symbol")
  } else if(method == "DESeq2") {
    table_DEGs <- listDEGs$DESeq2 magrittr::`%>%` tibble::rownames_to_column("Symbol")
  } else {
    table_DEGs <- listDEGs$edgeR magrittr::`%>%` tibble::rownames_to_column("Symbol")
  }

  if(padjusted) {
    table_DEGs <- table_DEGs magrittr::`%>%`
      tibble::rownames_to_column("Symbol") magrittr::`%>%`
      dplyr::mutate(Rank = log2FoldChange * -log10(padj)) magrittr::`%>%`
      dplyr::select(Symbol, Rank) magrittr::`%>%`
      dplyr::arrange(-Rank)
  } else {
    table_DEGs <- table_DEGs magrittr::`%>%`
      tibble::rownames_to_column("Symbol") magrittr::`%>%`
      dplyr::mutate(Rank = log2FoldChange * -log10(pvalue)) magrittr::`%>%`
      dplyr::select(Symbol, Rank) magrittr::`%>%`
      dplyr::arrange(-Rank)
  }

  return(table_DEGs)
}
