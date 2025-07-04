### Run infiltR ----------------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for infiltR, a tool for tumor microenvironmnet characterization
#' @note 2024-11-13
#' @title infiltR
#' @details
#'
#' InfiltR characterize the tumor microenvironment by evaluating the amount of
#' infiltrating immune cells of RNA-seq samples of solid tumor biopsies.
#' It relies on three complentary tools:
#' - [MCP-counter](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5);
#' - [CIBERSORT](https://pmc.ncbi.nlm.nih.gov/articles/PMC5895181/);
#' - [quanTIseq](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0638-6).
#'
#' InfiltR provides absolute and relative estimates of infiltrating immune cells,
#' it returns abundance tables for each tool that can be used for further downstream
#' analyses, and prints publication ready figures.
#'
#'
#' @param counts_table A RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param metadata A metadata matrix with samples as rows
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
#' @param sample_levels Order of samples groups for plotting.
#'   Default: "default", groups are plotted in alphabetical order.
#' @param my_palette A vector of color to be passed to the ggplot functions. If
#'   provided by the user, it must be the same length of number of factors in
#'   sample_groups.
#'   Default: "default", it uses standard ggplot2 palette.
#' @param plot_stats Boolean to control plotting of p-values brackets above MCP-
#'   counter absolute quantification violin plots.
#'   Default: TRUE
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
#' @param mcp_featuresType MCP-counter parameter.
#'   featuresType of the counts_table to be passed to MCP-counter functions.
#'   Possible values: "affy133P2_probesets", "HUGO_symbols", "ENTREZ_ID", "ENSEMBL_ID".
#'   Default: "HUGO_symbols"
#' @param mcp_probesets MCP-counter parameter.
#'   The probe set that indicates the genes of the immune infiltrating cells
#'   Default: "default", which retrieves the original MCP-counter repo with
#'   \code{read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character")}
#' @param mcp_genes MCP-counter parameter.
#'   The gene set that indicates the genes of the immune infiltrating cells
#'   Default: "default", which retrieves the original MCP-counter repo with
#'   \code{read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE}
#' @param cb_perm CIBERSORT parameter.
#'   Number of permutations to be performed to get a p-value on the estimated
#'   infiltrating immune cells.
#'   Default: 100
#' @param cb_pval_thr CIBERSORT parameter.
#'   P-value threshold to retain samples for comparisons.
#'   Default: 0.1
#' @param cb_QN CIBERSORT parameter.
#'   Perform Quantile Normalization on the counts_table.
#'   Default: FALSE
#' @param qs_is_arraydata quanTIseq parameter
#'   If it is microarray or RNAseq data.
#'   Default: FALSE
#' @param qs_method quanTIseq parameter
#'   If it is solid tumor data or PBMC data.
#'   Default: TRUE
#' @param qs_scale_mRNA quanTIseq parameter
#'   Logical value. If set to FALSE, it disables the correction of cell-type-specific
#'   mRNA content bias.
#'   Default: TRUE
#' @param qs_method quanTIseq parameter
#'   Character string, defining the deconvolution method to be used: lsei for
#'   constrained least squares regression, hampel, huber, or bisquare for robust
#'   regression with Huber, Hampel, or Tukey bisquare estimators, respectively.
#'   Default: "lsei".
#' @param qs_rm_genes quanTIseq parameter
#'   Character vector, specifying which genes have to be excluded from the
#'   deconvolution analysis. It can be provided as:
#'   • a vector of gene symbols (contained in the expression_data)
#'   • a single string among the choices of "none" (no genes are removed) and
#'     "default" (a list of genes with noisy expression RNA-seq data is removed,
#'     as explained in the quanTIseq paper).
#'   Default: "default" for RNA-seq data, "none" for microarrays.
#' @export
#'
#
# up_packages = c(
#   "BiocGenerics", "cowplot", "dplyr", "doSNOW", "e1071", "foreach", "ggplot2", "ggpmisc", "lubridate",
#   "MCPcounter", "parallel", "preprocessCore", "quantiseqr", "rlang", "rstatix", "scales", "stats", "stringr", "utils"
# )
# lapply(up_packages, require, character.only = TRUE)
infiltR = function(
    counts_table,
    metadata,
    sample_groups,
    sample_levels = "default",
    reference_samples = "default",
    my_palette = "default",
    plot_stats = TRUE,
    save_plots = TRUE,
    outdir = "default",
    mcp_featuresType = "HUGO_symbols",
    mcp_probesets = "default",
    mcp_genes = "default",
    cb_perm = 500,
    cb_pval_thr = 0.1,
    cb_QN = FALSE,
    qs_is_arraydata = FALSE,
    qs_is_tumordata = TRUE,
    qs_scale_mRNA = TRUE,
    qs_method = "lsei",
    qs_rm_genes = "default"
){

  # get starting time
  t0 = lubridate::now()
  message("[", as.POSIXct(lubridate::now()), "] ... Start infiltR!")

  # run estimates
  infiltr_out = run_estimates(
    counts_table,
    metadata,
    sample_groups,
    mcp_featuresType,
    mcp_probesets,
    mcp_genes,
    cb_perm = cb_perm,
    cb_pval_thr = cb_pval_thr,
    cb_QN = cb_QN,
    qs_is_arraydata = qs_is_arraydata,
    qs_is_tumordata = qs_is_tumordata,
    qs_scale_mRNA = qs_scale_mRNA,
    qs_method = qs_method,
    qs_rm_genes = qs_rm_genes
  )

  # print infiltration
  message("[", as.POSIXct(lubridate::now()), "] ... Generate figures")
  plot_infiltR(
    infiltr_out,
    metadata,
    sample_groups,
    sample_levels,
    my_palette,
    plot_stats,
    save_plots,
    outdir
  )

  # print QC report
  if(nrow(infiltr_out[["cb_sign"]]) >= 0){

    message("[", as.POSIXct(lubridate::now()), "] ... Generate CIBERSORT QC figures")
    plot_cb_qc(
      infiltr_out,
      metadata,
      sample_groups,
      save_plots,
      outdir
    )

  }

  # print estimates comparisons
  message("[", as.POSIXct(lubridate::now()), "] ... Generate comparisons figures")
  plot_comparisons(
    infiltr_out,
    save_plots,
    outdir
  )


  ### TODO: check if install
  ### TODO: add proper README and usage
  ### TODO: build package
  ### TODO: dependencies
  ### TODO: check function descriptions
  ### TODO: export function
  ### TODO: clean-up legacy code
  ### TODO: test / mock dataset


}
