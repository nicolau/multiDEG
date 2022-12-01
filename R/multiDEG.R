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

  # class.column <- "Class"
  # adjust.method <- "fdr"
  # treated <- "INF"
  # nontreated <- "CTRL"
  # covariables <- c("Gender")
  # covariables <- NULL
  # paired.samples.column = NULL
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
    phenodata <- phenodata %>% tidyr::unite(Concatenate, c(covariables), remove = FALSE)
    model <- paste0("~", paste(c(paired.samples.column, "Concatenate", class.column), collapse = "+"))
  } else if(length(covariables) == 1) {
    model <- paste0("~", paste(c(paired.samples.column, covariables, class.column), collapse = "+"))
  } else {
    model <- paste0("~", paste(c(paired.samples.column, class.column), collapse = "+"))
  }


  phenodata  <- phenodata %>% dplyr::filter(Class %in% c(treated, nontreated))
  raw.exp    <- raw.exp %>% dplyr::select(phenodata$Sample)

  conditions <- factor(t(phenodata$Class)) %>% relevel(conditions, ref = nontreated)

  # Use this variable as covariates
  # gender     <- factor(t(phenodata$Gender))

  if(!identical(colnames(raw.exp), phenodata$Sample)) {
    stop("Error: colnames(raw.exp) should be the same of phenodata$Sample")
  }

  ############################################################ Wilcoxon rank-sum test ############################################################
  message("Calculating DE genes using Wilcoxon rank-sum test...")
  # edgeR TMM normalize
  y <- DGEList(counts = raw.exp, group = conditions)
  ## Remove rows conssitently have zero or very low counts
  keep <- filterByExpr(y)
  y <- y[keep, keep.lib.sizes = FALSE]
  ## Perform TMM normalization and convert to CPM (Counts Per Million)
  y <- calcNormFactors(y, method = "TMM")
  norm.exp <- cpm(y) %>% as.data.frame()
  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(norm.exp), function(i){
    data <- cbind.data.frame(gene = as.numeric(t(norm.exp[i,])), conditions)
    p <- wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = adjust.method)
  # Calculate the fold-change for each gene
  dataCon1 <- norm.exp %>% dplyr::select(c(which(conditions==nontreated)))
  dataCon2 <- norm.exp %>% dplyr::select(c(which(conditions==treated)))
  foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  # Output results based on the FDR threshold 0.05
  outRst <- data.frame(row.names = rownames(norm.exp), log2foldChange = foldChanges, pValues = pvalues, FDR = fdr) %>% na.omit(outRst)

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

  tTags <- res %>% as.data.frame()
  result[['DESeq2']] <- tTags
  message("Done!")
  message("")
  #################################################################### DESeq2 ####################################################################




  ##################################################################### edgeR ####################################################################
  message("Calculating DE genes using edgeR...")
  # group  <- conditions
  y      <- DGEList(raw.exp)
  y      <- calcNormFactors(y)
  # term   <- ~0+group
  term   <- as.formula(model)
  design <- model.matrix(term)

  y      <- estimateDisp(y, design)
  fit    <- glmFit(y, design)
  lrt    <- glmLRT(fit, coef = 2, contrast = c(-1, 1))
  tTags  <- topTags(lrt, n = NULL)[["table"]]
  result[['edgeR']] <- tTags
  message("Done!")
  message("")
  ##################################################################### edgeR ####################################################################


  return(result)
}
