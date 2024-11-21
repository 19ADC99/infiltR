### Run MCP-counter ------------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for MCP-counter
#' @note 2024-11-13
#'
#' @param infiltr_out infiltR output object
#' @param counts_table A RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param metadata A metadata matrix with samples as rows
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
run_quantiseq = function(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    qs_is_arraydata,
    qs_is_tumordata,
    qs_scale_mRNA,
    qs_method,
    qs_rm_genes
){

  # get estimates
  message("[", as.POSIXct(lubridate::now()), "] ....... quanTIseq estimates")
  infiltr_out[["quantiseq"]] = quantiseqr::run_quantiseq(
    expression_data = counts_table,
    signature_matrix = "TIL10",
    is_arraydata = qs_is_arraydata,
    is_tumordata = qs_is_tumordata,
    scale_mRNA = qs_scale_mRNA,
    method = qs_method,
    rm_genes = qs_rm_genes
  )

  return(infiltr_out)

}
