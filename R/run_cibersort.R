### Run CIBERSORT --------------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for CIBERSORT
#' @note 2024-11-15
#'
#' @param infiltr_out infiltR output object
#' @param counts_table A RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param metadata A metadata matrix with samples as rows
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
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
run_cibersort <- function(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    cb_perm,
    cb_pval_thr,
    cb_QN
){

  # retrieve signature matrix
  message("[", as.POSIXct(lubridate::now()), "] ....... retrieve CIBERSORT signature matrix")
  if(mcp_probesets == "default"){
    mcp_probesets = read.table(
      curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),
      sep = "\t",
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
  }

  if(mcp_genes == "default"){
    mcp_genes = read.table(
      curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),
      sep = "\t",
      stringsAsFactors = FALSE,
      header = TRUE,
      colClasses = "character",
      check.names = FALSE
    )
  }

  # get estimates
  message("[", as.POSIXct(lubridate::now()), "] ....... MCP-counter estimates")
  infiltr_out[["mcp_counter"]] = MCPcounter::MCPcounter.estimate(
    counts_table,
    featuresType = "HUGO_symbols"
  ) %>%
    t() %>%
    as.data.frame()

  # total infiltrating
  infiltr_out[["mcp_counter"]] = infiltr_out[["mcp_counter"]] %>%
    dplyr::mutate(all_infiltrating = rowSums(.))
  infiltr_out[["mcp_counter"]]$group = metadata[[sample_groups]]

  return(infiltr_out)

}
