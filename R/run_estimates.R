### Run infiltR estimates ------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for infiltR estimates
#' @note 2024-12-11
#' @title run_estimates
#' @details
#'
#' This function runs CIBERSORT, MCP-counter and quanTIseq tools
#'
#' @param counts_table A RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param metadata A metadata matrix with samples as rows
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
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
run_estimates = function(
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
){

  # declare output obj
  infiltr_out = list()

  # run MCP-counter
  message("[", as.POSIXct(lubridate::now()), "] ... Run MCP-counter")
  infiltr_out = run_mcpcounter(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    mcp_featuresType,
    mcp_probesets,
    mcp_genes
  )

  # run CIBERSORT
  message("[", as.POSIXct(lubridate::now()), "] ... Run CIBERSORT")
  infiltr_out = run_cibersort(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    cb_perm = cb_perm,
    cb_pval_thr = cb_pval_thr,
    cb_QN = cb_QN
  )

  # run quanTIseq
  message("[", as.POSIXct(lubridate::now()), "] ... Run quanTIseq")
  infiltr_out = run_quantiseq(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    qs_is_arraydata = qs_is_arraydata,
    qs_is_tumordata = qs_is_tumordata,
    qs_scale_mRNA = qs_scale_mRNA,
    qs_method = qs_method,
    qs_rm_genes = qs_rm_genes
  )

  return(infiltr_out)

}
