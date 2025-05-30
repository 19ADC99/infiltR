% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cb_svm.R
\name{cb_svm}
\alias{cb_svm}
\title{cb_svm}
\usage{
cb_svm(sample, cb_LM22)
}
\arguments{
\item{sample}{A column of the RNAseq counts table, with genes as rows and samples as
columns. The count table is expected to be already preprocessed (i.e.: lowly
abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
normalized).}

\item{cb_LM22}{CIBERSORT parameter.
LM22 siganture matrix.

run CIBERSORT nu SVM}
}
\description{
cb_svm
}
\details{
Core Support Vector Machine algorithm of CIBERSORT.
}
\note{
run nu-SVM

2024-11-19
}
\author{
Andrea Del Cortona \href{mailto:andrea.delcortona@gmail.com}{andrea.delcortona@gmail.com}
}
