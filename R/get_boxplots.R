### Generate boxplots ----------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for p-values positions to be annotated
#' @note 2024-11-22
#' @title get_boxplots
#' @details
#'
#' This function generates boxplots for the main figure.
#'
#' @param long_matrix The long (melted) version of a cell adundance matrix, (e.g.:)
#'   as obtained from CIBERSORT or MCP-counter
#' @param ttest A t-test table obtained from rstatix::t_test
#' @param x feature in the x axis. Used as well to facet wrap if cell_type == TRUE
#' @param y feature in the y axis
#' @param fill feature to be used for the fill
#' @param my_y_positions for plotting p-values. My_y_positions vectors as
#'   obtained from the function get_pvalues_positions
#' @param my_x_min for plotting p-values. my_x_min vectors as
#'   obtained from the function get_pvalues_positions
#' @param my_x_max for plotting p-values. my_x_max vectors as
#'   obtained from the function get_pvalues_positions
#' @param my_annotation for plotting p-values. my_annotation vectors as
#'   obtained from the function get_pvalues_positions
#' @param my_palette A vector of color to be passed to the ggplot functions. If
#'   provided by the user, it must be the same length of number of factors in
#'   sample_groups.
#' @param plot_stats Boolean to control plotting of p-values brackets above MCP-
#'   counter absolute quantification violin plots.
#'   Default: TRUE
#' @param log_y y-axis to be plotted as log10. Useful to plot absolute counts (
#'   i.e.: as from MCP-counter).
#'   Default: FALSE
#' @param cell_type if TRUE, it expects a long_matrix with cell_types and it
#'   generates a faceted plot, once facet for each cell-type. Stylistic choice
#'   to make the cell type names readable
#' @export
get_boxplots = function(
    long_matrix,
    ttest,
    my_x,
    my_y,
    my_fill,
    my_y_positions,
    my_x_min,
    my_x_max,
    my_annotation,
    my_palette,
    plot_stats,
    log_y = FALSE,
    cell_type = TRUE
){

  # check if there were multiple comparisons
  adjusted = ifelse("p.adj.signif" %in% colnames(ttest), TRUE, FALSE)

  # general plot
  my_plot = ggplot(long_matrix, aes(x = .data[[my_x]], y = .data[[my_y]], fill = .data[[my_fill]])) +
    geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
    stat_summary(fun = "median", colour = "red", geom = "crossbar", size = 0.15,
                 position = position_dodge(0.9), width = 0.5) +
    scale_fill_manual(values = my_palette) +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank())

  # log scale?
  if(log_y == TRUE){
    my_plot = my_plot +
      scale_y_log10()
  }

  # add p-value if plot_stats == TRUE
  if(plot_stats){

    for(k in 1:length(my_annotation)){

      # check if there were multiple comparisons
      adjusted = ifelse("p.adj.signif" %in% colnames(ttest), TRUE, FALSE)

      my_plot = my_plot +
        ggplot2::annotate(
          "segment",
          x = my_x_min[[k]],
          xend = my_x_max[[k]],
          y = my_y_positions[[k]],
          yend = my_y_positions[[k]],
        )

      if(adjusted){

        my_plot = my_plot +
          ggplot2::annotate(
            "text",
            x = mean(c(my_x_min[[k]], my_x_max[[k]])),
            y = my_y_positions[[k]] + 0.0025,
            label = my_annotation[[k]]
          )

      } else {

        my_plot = my_plot +
          ggplot2::annotate(
            "text",
            x = mean(c(my_x_min[[k]], my_x_max[[k]])),
            y = my_y_positions[[k]]*1.5,
            label = my_annotation[[k]]
          )

      }

    }

  }


  # cell types or group?
  if(cell_type == TRUE){
    my_plot = my_plot +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title.y = element_blank(),
            legend.position = "none",
            panel.spacing = unit(0, "lines"))
  } else {
    my_plot = my_plot +
      ylab("cell counts") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "none")
  }

  return(my_plot)

}





