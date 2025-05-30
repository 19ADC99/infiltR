### Get p-values positions -----------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for p-values positions to be annotated
#' @note 2024-11-15
#' @title get_pvalues_positions
#' @details
#'
#' I get the x and y positions and the adj-p-value significance to be annotated
#' on ggplot boxplots for different cell types across the groups.
#' I replace similar functions in ggpubr and ggsign packages that seems not to
#' work well for me.
#'
#' @param long_matrix The long (melted) version of a cell adundance matrix, (e.g.:)
#'   as obtained from CIBERSORT or MCP-counter
#' @param ttest A t-test table obtained from rstatix::t_test
#' @param nfact number of distinct groups
#' @export
get_pvalues_positions = function(
    long_matrix,
    ttest,
    n_fact
){

  # get sign and positions
  my_y_positions = c()
  my_x_min = c()
  my_x_max = c()
  my_annotation = c()

  # check if there were multiple comparisons
  adjusted = ifelse("p.adj.signif" %in% colnames(ttest), TRUE, FALSE)

  # absolute VS cell-type p-values
  if("cell_type" %in% colnames(ttest)){

    for(k in 1:nrow(ttest)){

      if(k %% n_fact == 1){

        my_group1 = as.character(ttest[k, "group1"])
        my_group2 = as.character(ttest[k, "group2"])
        my_cell_type = as.character(ttest[k, "cell_type"][[1]])

        # check if absolute or relative cell abundance
        if(max(long_matrix$value) <= 1){

          my_y_val = long_matrix %>%
            dplyr::filter(cell_type == my_cell_type) %>%
            dplyr::select(value) %>%
            max()
          my_y_val = my_y_val + 0.025

        } else {

          my_y_val = long_matrix %>%
            dplyr::filter(cell_type == my_cell_type) %>%
            dplyr::select(value) %>%
            max() %>%
            round()
          my_y_val = my_y_val*1.4

        }

      } else {

        # check if absolute or relative cell abundance
        if(max(long_matrix$value) <= 1){
          my_y_val = my_y_positions[length(my_y_positions)] + 0.025
        } else {
          my_y_val = my_y_positions[length(my_y_positions)]*1.4
        }

      }

      my_y_positions = c(my_y_positions, my_y_val)
      my_x_min = c(my_x_min, as.numeric(ttest[k, "xmin"]))
      my_x_max = c(my_x_max, as.numeric(ttest[k, "xmax"]))

      # add significance level
      if(adjusted){

        if(as.character(ttest[k, "p.adj.signif"]) != "ns"){
          my_annotation = c(my_annotation, as.character(ttest[k, "p.adj.signif"]))
        } else {
          my_annotation = c(my_annotation, "")
        }

      } else {

        if(as.character(ttest[k, "p"]) != "ns"){
          my_annotation = c(my_annotation, as.character(ttest[k, "p"]))
        } else {
          my_annotation = c(my_annotation, "")
        }

      }

    }

  } else {

    for(k in 1:nrow(ttest)){

      if(k == 1){

        my_y_val = long_matrix %>%
          max() %>%
          round()
        my_y_val = my_y_val*1.4

      } else {

        my_y_val = my_y_positions[length(my_y_positions)]*1.4

      }

      my_y_positions = c(my_y_positions, my_y_val)
      my_x_min = c(my_x_min, as.numeric(ttest[k, "xmin"]))
      my_x_max = c(my_x_max, as.numeric(ttest[k, "xmax"]))

      # add significance level
      if(adjusted){

        if(as.character(ttest[k, "p.adj.signif"]) != "ns"){
          my_annotation = c(my_annotation, as.character(ttest[k, "p.adj.signif"]))
        } else {
          my_annotation = c(my_annotation, "")
        }

      } else {

        if(as.character(ttest[k, "p"]) != "ns"){
          my_annotation = c(my_annotation, as.character(ttest[k, "p"]))
        } else {
          my_annotation = c(my_annotation, "")
        }

      }

    }

  }


  return(
    list(
      "my_y_positions" = my_y_positions,
      "my_x_min" = my_x_min,
      "my_x_max" = my_x_max,
      "my_annotation" = my_annotation
    )
  )

}
