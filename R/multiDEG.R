#' Calculates DE genes table for Wilcoxon rank-sum test, DESeq2 and edgeR
#'
#' @param raw.exp data.frame of raw counts
#' @param phenodata data.frame of phenotypic data for all samples
#' @param treated A string containing the name of class to be considered as treatment
#' @param nontreated A string containing the name of class to be considered as nontreatment
#' @param class.column A string contatining the name of column where treated and nontreated will be extracted
#' @param adjust.method Merhod of multiple comparison will be used for Wilcoxon rank-sum test pvalues.
#' @param covariables A vector list of variables should be used to covariates adjusment.
#' @param paired.samples.column A boolean variable to indicate paired-sample comparison
#'
#' @return A list of data.frame for Wilcoxon rank-sum test, DESeq2 and edgeR
#' @export
DEG_analysis <- function(raw.exp, phenodata, treated, nontreated, class.column = "Class",
                         adjust.method = c( "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                         covariables = NULL, paired.samples.column = NULL) {

  # class.column          <- "Class"
  # adjust.method         <- "fdr"
  # treated               <- "INF"
  # nontreated            <- "CTRL"
  # covariables           <- c("Gender")
  # covariables           <- NULL
  # paired.samples.column <- NULL
  # data(phenodata)
  # data(raw.exp)

  if(!require(tidyverse)) { stop("tidyverse package not available.") }
  if(!require(DESeq2)) { stop("DESeq2 package not available.") }
  if(!require(edgeR)) { stop("edgeR package not available.") }
  if(!require(rstatix)) { stop("rstatix package not available.") }
  if(!require(grid)) { stop("rstatix package not available.") }
  if(!require(ashr)) { stop("ashr package not available.")}


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

  if(!class.column %in% colnames(phenodata)) {stop("Class column: \"", paste(class.column, collapse = ", "), "\" has not found in the phenodata file!")}

  if(class.column != "Class") {
    if(sum(colnames(phenodata) == "Class") == 1) {
      colnames(phenodata)[colnames(phenodata) == "Class"] <- "ClassOld"
    }
    colnames(phenodata)[colnames(phenodata) == class.column] <- "Class"
  }


  if(length(covariables) > 1) {
    phenodata <- phenodata %>% tidyr::unite(Concatenate, c(covariables), remove = FALSE)
    model <- paste0("~0+", paste(c(paired.samples.column, "Concatenate", class.column), collapse = "+"))
  } else if(length(covariables) == 1) {
    model <- paste0("~0+", paste(c(paired.samples.column, covariables, class.column), collapse = "+"))
  } else {
    model <- paste0("~0+", paste(c(paired.samples.column, class.column), collapse = "+"))
  }


  phenodata  <- phenodata %>% dplyr::filter(Class %in% c(treated, nontreated))
  raw.exp    <- raw.exp %>% dplyr::select(phenodata$Sample)

  conditions <- factor(t(phenodata$Class)) %>% relevel(conditions, ref = nontreated)

  if(!identical(colnames(raw.exp), phenodata$Sample)) { stop("Error: colnames(raw.exp) should be the same of phenodata$Sample") }






  ############################################################ Wilcoxon rank-sum test ############################################################
  message("Calculating DE genes using Wilcoxon rank-sum test...")
  # edgeR TMM normalize
  y <- edgeR::DGEList(counts = raw.exp, group = conditions)
  ## Remove rows conssitently have zero or very low counts
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, keep.lib.sizes = FALSE]
  ## Perform TMM normalization and convert to CPM (Counts Per Million)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  norm.exp <- edgeR::cpm(y) %>% as.data.frame()


  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(norm.exp), function(i){
    data <- cbind.data.frame(gene = as.numeric(t(norm.exp[i,])), conditions) #, gender = phenodata$Gender)
    p <- wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = adjust.method)
  # Calculate the fold-change for each gene
  dataCon1 <- norm.exp %>% dplyr::select(c(which(conditions==nontreated)))
  dataCon2 <- norm.exp %>% dplyr::select(c(which(conditions==treated)))
  foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  # Output results based on the FDR threshold 0.05
  outRst <- data.frame(row.names = rownames(norm.exp), log2FoldChange = foldChanges, pvalue = pvalues, padj = fdr) %>% na.omit(outRst)

  result[['Wilcoxon']] <- outRst
  message("Done!")
  message("")
  ############################################################ Wilcoxon rank-sum test ############################################################






  #################################################################### limma #####################################################################
  message("Calculating DE genes using limma...")
  design <- model.matrix(as.formula(model))
  fit <- lmFit(log2(edgeR::cpm(raw.exp)+1), design)
  contrast_formula <- paste0("Class", treated, " - Class", nontreated, collapse = "")
  contrasts <- makeContrasts(contrast_formula, levels = design)

  fit2 <- contrasts.fit(fit = fit, contrasts = contrasts)
  fit2 <- eBayes(fit2)
  tTags <- topTable(fit = fit2, coef = contrast_formula, number = Inf, adjust.method = "BH") %>%
    as.data.frame() %>%
    dplyr::rename(log2FoldChange = logFC,
                  baseMean       = AveExpr,
                  pvalue         = P.Value,
                  padj           = adj.P.Val)

  result[['limma']] <- tTags
  message("Done!")
  message("")
  #################################################################### limma #####################################################################






  #################################################################### DESeq2 ####################################################################
  message("Calculating DE genes using DESeq2...")
  contrasts <- c("Class", treated, nontreated)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw.exp, colData = phenodata, design = as.formula(model))

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrasts, name = paste0(contrasts[2], "_vs_", contrasts[3]))
  # res_tableOE <- lfcShrink(dds, coef = resultsNames(dds), type = "apeglm")

  # resultsNames(dds) %in% resultsNames(dds)
  # res_tableOE <- lfcShrink(dds, coef = resultsNames(dds)[1], res = res, type = "apeglm")
  # res_tableOE <- lfcShrink(dds, contrast = contrasts, res = res, type = "ashr")
  # res_tableOE <- lfcShrink(dds, contrast = contrasts, res = res, type = "normal")

  # tTags <- res_tableOE %>% as.data.frame()
  #
  # ggplot(data = tTags, aes(x = log10(baseMean), y = log2FoldChange)) +
  #   geom_hline(yintercept = 1) +
  #   geom_hline(yintercept = -1) +
  #   geom_point()
  #
  # ggplot(data = tTags, aes(x = log10(baseMean), y = log2FoldChange)) +
  #   geom_point()

  tTags <- res %>% as.data.frame()
  result[['DESeq2']] <- tTags
  message("Done!")
  message("")
  #################################################################### DESeq2 ####################################################################






  ##################################################################### edgeR ####################################################################
  message("Calculating DE genes using edgeR...")
  if(length(covariables) > 0) {
    phenodata <- phenodata %>% tidyr::unite(Concatenate, c(covariables), remove = FALSE)
    model <- paste0("~0+", paste(c(paired.samples.column, paste0("Cov", 1:length(covariables)), class.column), collapse = "+"))
    contrast_vector <- c(-1, 1, 0)
    for(i in length(covariables)) {
      assign(paste0("Cov", i), phenodata[, covariables[i]] %>% unlist(use.names = F))
    }
  }else if(!is.null(paired.samples.column)) {
    model <- paste0("~0+", paste(c(paired.samples.column, class.column), collapse = "+"))
    contrast_vector <- c(-1, 1, 0)
  } else {
    model <- paste0("~0+", paste(c(paired.samples.column, class.column), collapse = "+"))
    contrast_vector <- c(-1, 1)
  }

  Class  <- conditions
  y      <- edgeR::DGEList(raw.exp)
  y      <- edgeR::calcNormFactors(y)
  term   <- as.formula(model)
  design <- model.matrix(term)

  y      <- edgeR::estimateDisp(y, design)
  fit    <- edgeR::glmFit(y, design)
  lrt    <- edgeR::glmLRT(fit, coef = 2, contrast = contrast_vector)
  tTags  <- edgeR::topTags(lrt, n = NULL)[["table"]]
  result[['edgeR']] <- tTags
  message("Done!")
  message("")
  ##################################################################### edgeR ####################################################################

  result[['edgeR']]   <- result[['edgeR']] %>% dplyr::rename(log2FoldChange = logFC, pvalue = PValue, padj = FDR)

  ################################################################# Overlap DEGs #################################################################
  overlap_result      <- overlap_DEGs(listDEGs = result)
  result[["overlap"]] <- overlap_result
  ################################################################# Overlap DEGs #################################################################

  return(result)
}


#' Calculates overlapping genes between Wilcox rank-sum test, DESeq2 and edgeR DE methods
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of Venn Diagram plots for up- and down-regulated genes
#' @export
overlap_DEGs_old <- function(listDEGs, p_cutoff = 0.05, log2fc_cutoff = 1, padjusted = F) {

  genes_wilcox_up   <- listDEGs$Wilcoxon %>% tibble::rownames_to_column("Symbol")
  genes_wilcox_down <- listDEGs$Wilcoxon %>% tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_wilcox_up <- genes_wilcox_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)

    genes_wilcox_down <- genes_wilcox_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  } else {
    genes_wilcox_up <- genes_wilcox_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)

    genes_wilcox_down <- genes_wilcox_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  }


  genes_edgeR_up   <- listDEGs$edgeR %>% tibble::rownames_to_column("Symbol")
  genes_edgeR_down <- listDEGs$edgeR %>% tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_edgeR_up<- genes_edgeR_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
    genes_edgeR_down <- genes_edgeR_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  } else {
    genes_edgeR_up<- genes_edgeR_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
    genes_edgeR_down <- genes_edgeR_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  }


  genes_DESeq2_up   <- listDEGs$DESeq2 %>% tibble::rownames_to_column("Symbol")
  genes_DESeq2_down <- listDEGs$DESeq2 %>% tibble::rownames_to_column("Symbol")
  if(padjusted) {
    genes_DESeq2_up <- genes_DESeq2_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
    genes_DESeq2_down <- genes_DESeq2_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  } else {
    genes_DESeq2_up <- genes_DESeq2_up %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff  & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
    genes_DESeq2_down <- genes_DESeq2_down %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff & pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  }


  plot_down <- plot.triple.venn(a1 = genes_wilcox_down, a2 = genes_edgeR_down, a3 = genes_DESeq2_down, labels = c("Wilcox", "edgeR", "DESeq2"))
  plot_up   <- plot.triple.venn(a1 = genes_wilcox_up,   a2 = genes_edgeR_up,   a3 = genes_DESeq2_up,   labels = c("Wilcox", "edgeR", "DESeq2"))

  if(!is.null(dev.list())) { dev.off() }

  return(list("plot_up" = plot_up, "plot_down" = plot_down))

}


#' Calculates overlapping genes between Wilcox rank-sum test, DESeq2 and edgeR DE methods
#'
#' @param listDEGs list of results from DEG_analysis
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A list of Venn Diagram plots for up- and down-regulated genes
#' @export
overlap_DEGs <- function(listDEGs, p_cutoff = 0.05, log2fc_cutoff = 1, padjusted = F) {

  degs <- list()

  for(method in names(listDEGs)) {
    genes_up   <- listDEGs[[method]] %>% tibble::rownames_to_column("Symbol")
    genes_down <- listDEGs[[method]] %>% tibble::rownames_to_column("Symbol")

    if(padjusted) {
      genes_up <- genes_up %>%
        dplyr::filter(log2FoldChange > log2fc_cutoff  & padj < p_cutoff) %>%
        dplyr::select(Symbol) %>%
        unlist(use.names = F)

      genes_down <- genes_down %>%
        dplyr::filter(log2FoldChange < -log2fc_cutoff & padj < p_cutoff) %>%
        dplyr::select(Symbol) %>%
        unlist(use.names = F)
    } else {
      genes_up <- genes_up %>%
        dplyr::filter(log2FoldChange > log2fc_cutoff  & pvalue < p_cutoff) %>%
        dplyr::select(Symbol) %>%
        unlist(use.names = F)

      genes_down <- genes_down %>%
        dplyr::filter(log2FoldChange < -log2fc_cutoff & pvalue < p_cutoff) %>%
        dplyr::select(Symbol) %>%
        unlist(use.names = F)
    }
    degs[[method]][["up"]]   <- genes_up
    degs[[method]][["down"]] <- genes_down

  }


  plot_down <- plot.triple.venn(
    a1 = degs[["limma"]][["down"]],
    a2 = degs[["edgeR"]][["down"]],
    a3 = degs[["DESeq2"]][["down"]],
    labels = c("limma", "edgeR", "DESeq2"))
  plot_up   <- plot.triple.venn(
    a1 = degs[["limma"]][["up"]],
    a2 = degs[["edgeR"]][["up"]],
    a3 = degs[["DESeq2"]][["up"]],
    labels = c("limma", "edgeR", "DESeq2"))

  if(!is.null(dev.list())) { dev.off() }

  return(list("plot_up" = plot_up, "plot_down" = plot_down))

}


#' Plot triple venn diagram
#'
#' @param a1 elements for the group 1
#' @param a2 elements for the group 2
#' @param a3 elements for the group 3
#' @param labels labels for each groups
#'
#' @return A list contains a Venn Diagram plot, and list of genes overlapped and exclusive for each group
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
#' @param type plot up- or down-regulated genes
#'
#' @export
plot <- function(listDEGs, type = c("up", "down")) {
  if(!is.null(dev.list())) { dev.off() }
  if(type == "up") {
    grid::grid.draw(listDEGs$overlap$plot_up$plot)
  } else {
    grid::grid.draw(listDEGs$overlap$plot_down$plot)
  }
}


#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param method DE method
#'
#' @return A complete list of genes, log2FC, pvalue, and padj
#' @export
get_DEG_table <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR")) {
  table_DEGs <- NULL
  if(method == "Wilcox") {
    table_DEGs <- listDEGs$Wilcoxon
  } else {
    table_DEGs <- listDEGs[[method]]
  }
  return(table_DEGs)
}


#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param method DE method
#' @param direction Direction of differential expression (up or down)
#' @param p_cutoff cut-off for p-value
#' @param log2fc_cutoff cut-off for log2 Fold-change values
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A vector of symbols
#' @export
get_DEG_symbols <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), direction = c("up", "down"),  p_cutoff = 0.05, log2fc_cutoff = 1, padjusted = F) {

  symbols    <- NULL
  table_DEGs <- NULL

  if(method == "Wilcox") {
    table_DEGs <- listDEGs$Wilcoxon
  } else if(method == "DESeq2") {
    table_DEGs <- listDEGs$DESeq2
  } else {
    table_DEGs <- listDEGs$edgeR
  }

  if(direction == "up") {
    symbols <- table_DEGs %>%
      dplyr::filter(log2FoldChange > log2fc_cutoff) %>%
      tibble::rownames_to_column("Symbol")
  } else {
    symbols <- table_DEGs %>%
      dplyr::filter(log2FoldChange < -log2fc_cutoff) %>%
      tibble::rownames_to_column("Symbol")
  }

  if(padjusted) {
    symbols <- symbols %>%
      dplyr::filter(padj < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  } else {
    symbols <- symbols %>%
      dplyr::filter(pvalue < p_cutoff) %>%
      dplyr::select(Symbol) %>%
      unlist(use.names = F)
  }

  return(symbols)
}

#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param method DE method
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A histogram plot for pvalues or padj
#' @export
pvalue_distribution_plot <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), padjusted = F) {

  table_DEGs <- NULL

  if(method == "wilcox") {
    table_DEGs <- listDEGs$Wilcoxon
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
#' @param method_one DE method for x-axis
#' @param method_two DE method for y-axis
#'
#' @return A scatter plot comparing two different DE methods
#' @export
log2fc_correlation_plot <- function(listDEGs, method_one = c("Wilcox", "DESeq2", "edgeR"), method_two = c("Wilcox", "DESeq2", "edgeR")) {

  table_one <- NULL
  table_two <- NULL

  if(method_one != method_two) {
    if(method_one == "Wilcox") {
      table_one <- listDEGs$Wilcoxon %>% tibble::rownames_to_column("Symbol")
      xlabel <- "Wilcox"
    } else if(method_one == "DESeq2") {
      table_one <- listDEGs$DESeq2 %>% tibble::rownames_to_column("Symbol")
      xlabel <- "DESeq2"
    } else {
      table_one <- listDEGs$edgeR %>% tibble::rownames_to_column("Symbol")
      xlabel <- "edgeR"
    }

    if(method_two == "Wilcox") {
      table_two <- listDEGs$Wilcoxon %>% tibble::rownames_to_column("Symbol")
      ylabel <- "Wilcox"
    } else if(method_two == "DESeq2") {
      table_two <- listDEGs$DESeq2 %>% tibble::rownames_to_column("Symbol")
      ylabel <- "DESeq2"
    } else {
      table_two <- listDEGs$edgeR %>% tibble::rownames_to_column("Symbol")
      ylabel <- "edgeR"
    }

    table_merge <- dplyr::inner_join(x = table_one, y = table_two, by = "Symbol")

    p <- ggplot2::ggplot(table_merge, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
      ggplot2::geom_point() +
      ggplot2::theme_classic() +
      xlab(xlabel) +
      ylab(ylabel)

    return(p)

  } else {
    stop("Please, select different methods for method_one and method_two parameters.")
  }
}


#' Export DEG tables for a specific DE method
#'
#' @param listDEGs list of results from DEG_analysis
#' @param method DE method
#' @param padjusted Boolean value to use or not adjusted p-values for p_cutoff
#'
#' @return A ranked table for GSEA tools
#' @export
save_ranked_table <- function(listDEGs, method = c("Wilcox", "DESeq2", "edgeR"), padjusted = FALSE) {

  table_DEGs <- NULL

  if(method == "wilcox") {
    table_DEGs <- listDEGs$Wilcoxon %>% tibble::rownames_to_column("Symbol")
  } else if(method == "DESeq2") {
    table_DEGs <- listDEGs$DESeq2 %>% tibble::rownames_to_column("Symbol")
  } else {
    table_DEGs <- listDEGs$edgeR %>% tibble::rownames_to_column("Symbol")
  }

  if(padjusted) {
    table_DEGs <- table_DEGs %>%
      # tibble::rownames_to_column("Symbol") %>%
      dplyr::mutate(Rank = log2FoldChange * -log10(padj)) %>%
      dplyr::select(Symbol, Rank) %>%
      dplyr::arrange(-Rank)
  } else {
    table_DEGs <- table_DEGs %>%
      # tibble::rownames_to_column("Symbol") %>%
      dplyr::mutate(Rank = log2FoldChange * -log10(pvalue)) %>%
      dplyr::select(Symbol, Rank) %>%
      dplyr::arrange(-Rank)
  }

  return(table_DEGs)
}
