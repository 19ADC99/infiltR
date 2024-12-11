### CIBERSORT permutations -----------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note run N permutations to get the null distribution
#' @note 2024-11-19
#'
#' @param counts_table A RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param cb_LM22 CIBERSORT parameter.
#'   LM22 siganture matrix.
#' @param cb_perm CIBERSORT parameter.
#'   Number of permutations to be performed to get a p-value on the estimated
#'   infiltrating immune cells.
#'   Default: 500
#'
#' run N (cb_perm) permutations to get the null distribution
get_cb_perm = function(
    counts_table,
    cb_LM22,
    cb_perm
    ){

  # declare output
  gene_expr_values = BiocGenerics::unlist(counts_table)
  distance = c()

  # run permutations
  pb = utils::txtProgressBar(min = 0, max = cb_perm, style = 3, file = stderr())
  for(k in 1:cb_perm){

    utils::setTxtProgressBar(pb, k)

    # get random expression
    random_expr = sample(gene_expr_values, nrow(cb_LM22))
    random_expr_scaled = (random_expr - mean(random_expr)) / stats::sd(random_expr)
    cb_estimates = cb_svm(random_expr_scaled, cb_LM22)
    distance[k] = as.numeric(cb_estimates$bm_r)

  }
  close(pb)

  return(distance)

}


