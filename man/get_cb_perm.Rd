% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cb_perm.R
\name{get_cb_perm}
\alias{get_cb_perm}
\title{get_cb_perm}
\usage{
get_cb_perm(counts_table, cb_LM22, cb_perm)
}
\arguments{
\item{counts_table}{A RNAseq counts table, with genes as rows and samples as
columns. The count table is expected to be already preprocessed (i.e.: lowly
abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
normalized).}

\item{cb_LM22}{CIBERSORT parameter.
LM22 siganture matrix.}

\item{cb_perm}{CIBERSORT parameter.
Number of permutations to be performed to get a p-value on the estimated
infiltrating immune cells.
Default: 500}
}
\description{
get_cb_perm
}
\details{
This function runs N (cb_perm) permutations to get the null distribution.
}
\note{
run N permutations to get the null distribution

2024-11-19
}
\author{
Andrea Del Cortona \href{mailto:andrea.delcortona@gmail.com}{andrea.delcortona@gmail.com}
}
