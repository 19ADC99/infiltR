### Run CIBERSORT --------------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for CIBERSORT
#' @note 2024-11-15
#'
#' Retro-engineered from https://www.nature.com/articles/nmeth.3337
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
#'   Default: 500
#' @param cb_pval_thr CIBERSORT parameter.
#'   P-value threshold to retain samples for comparisons.
#'   Default: 0.1
#' @param cb_QN CIBERSORT parameter.
#'   Perform Quantile Normalization on the counts_table.
#'   Default: FALSE
run_cibersort = function(
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
  cb_LM22 = read.table(
    curl:::curl("https://raw.githubusercontent.com/19ADC99/infiltR/master/data/cb_LM22.txt"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  # counts table quantile normalization
  message("[", as.POSIXct(lubridate::now()), "] ....... run quantile normalization")
  if(cb_QN == TRUE){

    sample_names = colnames(counts_table)
    gene_names = rownames(counts_table)

    counts_table = counts_table %>%
      as.matrix() %>%
      preprocessCore::normalize.quantiles() %>%
      BiocGenerics::as.data.frame()

    colnames(counts_table) = sample_names
    rownames(counts_table) = gene_names

  }

  # limit the analysis to genes in LM22 for increased speed
  counts_table_red = counts_table %>%
    dplyr::filter(row.names(counts_table) %in% row.names(cb_LM22))
  cb_LM22 = cb_LM22 %>%
    dplyr::filter(row.names(cb_LM22) %in% row.names(counts_table_red))


  # standardize LM22 matrix
  cb_LM22 = (as.matrix(cb_LM22) - mean(as.matrix(cb_LM22))) / stats::sd(as.matrix(cb_LM22))

  # empirical null distribution of correlation coefficients
  message("[", as.POSIXct(lubridate::now()), "] ....... estimate empirical null distribution with permutations")
  if(cb_perm > 0){

    nulldist = get_cb_perm(counts_table_red, cb_LM22, cb_perm) %>% sort()

  }

  # create output table
  infiltr_out[["cb_all"]] = matrix(nrow = 0, ncol = 25) %>%
    data.frame() %>%
    stats::setNames(c(colnames(cb_LM22), "P-value", "Correlation", "RMSE"))

  # iterate samples
  message("[", as.POSIXct(lubridate::now()), "] ....... deconvolute samples")
  pb = utils::txtProgressBar(min = 0, max = ncol(counts_table_red), style = 3, file = stderr())
  for(k in 1:ncol(counts_table_red)){

    utils::setTxtProgressBar(pb, k)

    # normalize sample
    sample = counts_table_red[, k]
    #sample = (sample - mean(sample)) / stats::sd(sample)

    # get estimates
    cb_results = cb_svm(sample, cb_LM22)

    # calculate p-value against empirical null distribution
    if (cb_perm > 0) {
      pvalue = 1 - (which.min(abs(nulldist - cb_results$bm_r)) / length(nulldist))
    } else {
      pvalue = 1
    }

    # append sample row
    infiltr_out[["cb_all"]][k, ] = c(
      cb_results$w,
      pvalue,
      cb_results$bm_rmse,
      cb_results$bm_r
    )
    rownames(infiltr_out[["cb_all"]])[k] = colnames(counts_table)[k]

  }
  close(pb)


  # cb table only significant samples
  infiltr_out[["cb_sign"]] = infiltr_out[["cb_all"]] %>%
    dplyr::filter(`P-value` <= cb_pval_thr)

  return(infiltr_out)

}
