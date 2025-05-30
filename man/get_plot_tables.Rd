% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_plot_tables.R
\name{get_plot_tables}
\alias{get_plot_tables}
\title{get_plot_tables}
\usage{
get_plot_tables(infiltr_out)
}
\arguments{
\item{infiltr_out}{infiltR output object}

\item{metadata}{A metadata matrix with samples as rows}

\item{sample_groups}{The columns of the metadata matrix to be used for the
comparisons between sample groups. E.g.: a column indicating which sample
was treated and which one not.}

\item{my_palette}{A vector of color to be passed to the ggplot functions. If
provided by the user, it must be the same length of number of factors in
sample_groups.
Default: "default", it uses standard ggplot2 palette.}

\item{save_plots}{save plots in pdf and png format
Default: TRUE}

\item{outdir}{Output directory of the plots
Default: "default", it prints in the current working directory}
}
\description{
get_plot_tables
}
\details{
This function prepr estimates tables for plotting
}
\note{
wrapper for plotting infiltR summary figures

2024-11-13
}
\author{
Andrea Del Cortona \href{mailto:andrea.delcortona@gmail.com}{andrea.delcortona@gmail.com}
}
