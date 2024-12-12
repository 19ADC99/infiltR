### CIBERSORT algorithm --------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note run nu-SVM
#' @note 2024-11-19
#' @title cb_svm
#' @details
#'
#' Core Support Vector Machine algorithm of CIBERSORT.
#'
#' @param sample A column of the RNAseq counts table, with genes as rows and samples as
#'   columns. The count table is expected to be already preprocessed (i.e.: lowly
#'   abundant genes removed with \code{edgeR::filterByExpr()} and \code{cpm()}
#'   normalized).
#' @param cb_LM22 CIBERSORT parameter.
#'   LM22 siganture matrix.
#'
#' run CIBERSORT nu SVM
#' @export
cb_svm = function(
    sample,
    cb_LM22
    ){

  # create CPU cluster
  cl = parallel::makeCluster(parallel::detectCores()[1] - 1, type = "SOCK")
  doSNOW::registerDoSNOW(cl)

  # run SVM with 3 different nu values
  model_output = foreach::foreach(k = 1:3) %dopar% {

    # get nu
    if(k == 1){ my_nu = 0.25 }
    if(k == 2){ my_nu = 0.5 }
    if(k == 3){ my_nu = 0.75 }

    # run SVM
    model = e1071::svm(
      cb_LM22,
      sample,
      type = "nu-regression",
      kernel = "linear",
      nu = my_nu,
      scale = FALSE
    )

    model

  }

  parallel::stopCluster(cl)


  # calculate weights
  rmses = c()
  corrv = c()
  for(k in 1:3){

    weights = t(model_output[[k]]$coefs) %*% model_output[[k]]$SV
    weights[which(weights < 0)] = 0
    rel_weigth = weights/sum(weights)
    multiplied_w = sweep(cb_LM22, MARGIN = 2, rel_weigth, "*")
    sum_w <- apply(multiplied_w, 1, sum)
    rmses[k] = sqrt(mean((sum_w - sample)^2))
    corrv[k] = stats::cor(sum_w, sample)

  }

  # select best model
  best = BiocGenerics::which.min(rmses)
  best_model = model_output[[best]]

  # get and normalize coefficients
  bm_coeff = t(best_model$coefs) %*% best_model$SV
  bm_coeff[which(bm_coeff < 0)] = 0
  w = (bm_coeff/sum(bm_coeff))
  bm_rmse = rmses[best]
  bm_r = corrv[best]

  cb_results = list("w" = w, "bm_rmse" = bm_rmse, "bm_r" = bm_r)

  return(cb_results)

}
