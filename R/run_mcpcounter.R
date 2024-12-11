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
run_mcpcounter = function(
    infiltr_out,
    counts_table,
    metadata,
    sample_groups,
    mcp_featuresType,
    mcp_probesets,
    mcp_genes
){

  # retrieve signature matrix
  message("[", as.POSIXct(lubridate::now()), "] ....... retrieve MCP-counter signature matrix")
  if(length(mcp_probesets) == 1){
    if(mcp_probesets == "default"){
      mcp_probesets = read.table(
        curl:::curl("https://raw.githubusercontent.com/19ADC99/infiltR/master/data/mcp_probesets.txt"),
        sep = "\t",
        stringsAsFactors = FALSE,
        colClasses = "character"
      )
    }
  }

  if(length(mcp_genes) == 1){
    if(mcp_genes == "default"){
      mcp_genes = read.table(
        curl:::curl("https://raw.githubusercontent.com/19ADC99/infiltR/master/data/mcp_genes.txt"),
        sep = "\t",
        stringsAsFactors = FALSE,
        header = TRUE,
        colClasses = "character",
        check.names = FALSE
      )
    }
  }

  # get estimates
  message("[", as.POSIXct(lubridate::now()), "] ....... MCP-counter estimates")
  infiltr_out[["mcp_counter"]] = MCPcounter::MCPcounter.estimate(
    counts_table,
    featuresType = "HUGO_symbols"
  ) %>%
    t() %>%
    as.data.frame()


  # normalized infiltrating
  infiltr_out[["mcp_counter_norm"]] = infiltr_out[["mcp_counter"]] %>%
    dplyr::mutate(all_infiltrating = rowSums(.)) %>%
    dplyr::mutate_all(~ ./all_infiltrating) %>%
    dplyr::select(-all_infiltrating)
  infiltr_out[["mcp_counter_norm"]]$group = metadata[[sample_groups]]

  # total infiltrating
  infiltr_out[["mcp_counter"]] = infiltr_out[["mcp_counter"]] %>%
    dplyr::mutate(all_infiltrating = rowSums(.))
  infiltr_out[["mcp_counter"]]$group = metadata[[sample_groups]]

  return(infiltr_out)

}
