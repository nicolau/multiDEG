% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiDEG.R
\name{get_DEG_symbols}
\alias{get_DEG_symbols}
\title{Export DEG tables for a specific DE method}
\usage{
get_DEG_symbols(
  listDEGs,
  method = c("Wilcox", "DESeq2", "edgeR"),
  direction = c("up", "down"),
  p_cutoff = 0.05,
  log2fc_cutoff = 1,
  padjusted = F
)
}
\arguments{
\item{listDEGs}{list of results from DEG_analysis}

\item{method}{DE method}

\item{direction}{Direction of differential expression (up or down)}

\item{p_cutoff}{cut-off for p-value}

\item{log2fc_cutoff}{cut-off for log2 Fold-change values}

\item{padjusted}{Boolean value to use or not adjusted p-values for p_cutoff}
}
\value{
A vector of symbols
}
\description{
Export DEG tables for a specific DE method
}
