### Plot infiltR figures -------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for plotting infiltR summary figures
#' @note 2024-11-13
#' @title plot_infiltR
#' @details
#'
#' This function plots infiltR estimates.
#'
#' @param infiltr_out infiltr_out output object
#' @param metadata A metadata matrix with samples as rows
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
#' @param my_palette A vector of color to be passed to the ggplot functions. If
#'   provided by the user, it must be the same length of number of factors in
#'   sample_groups.
#'   Default: "default", it uses standard ggplot2 palette.
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
#' @export
plot_infiltR = function(
    infiltr_out,
    metadata,
    sample_groups,
    my_palette,
    save_plots = TRUE,
    outdir = "default"
){

  #### General ####

  # get number of metadata groups
  n_fact = metadata[[sample_groups]] %>%
      as.factor() %>%
      levels() %>%
      length()

  # check my palette
  if(length(my_palette) == 1){
    if(my_palette == "default"){
      my_palette = scales::hue_pal()(n_fact)
    }
  }

  # prep dataframes for printing figures
  message("[", as.POSIXct(lubridate::now()), "] ... Prepping dataframes for figures and report")
  plot_tables = get_plot_tables(
    infiltr_out,
    metadata,
    sample_groups
  )


  #### MCP-counter ####

  mcp_absolute_infil = infiltr_out[["mcp_counter"]]

  # t-test
  tests_abs = mcp_absolute_infil %>%
    rstatix::t_test(all_infiltrating ~ group)
  tests_abs = tests_abs %>% rstatix::add_xy_position()

  my_pvalues = get_pvalues_positions(mcp_absolute_infil$all_infiltrating, tests_abs)

  p_mcp_abs_all = get_boxplots(
    mcp_absolute_infil,
    "group",
    "all_infiltrating",
    "group",
    my_y_positions = my_pvalues[["my_y_positions"]],
    my_x_min = my_pvalues[["my_x_min"]],
    my_x_max = my_pvalues[["my_x_max"]],
    my_annotation = my_pvalues[["my_annotation"]],
    my_palette,
    log_y = TRUE,
    cell_type = FALSE
  )

  p_mcp_abs_all = p_mcp_abs_all +
    labs(title = "Absolute quantification")


  # get transformed and plot per cell abundances
  mcp_absolute_infil = plot_tables[["mcp_absolute_infil"]]

  ## retrieve sign annotation and sign positions
  # t-test
  tests_abscell = mcp_absolute_infil %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_abscell = tests_abscell %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(mcp_absolute_infil, tests_abscell)

  p_mcp_abs_cells = get_boxplots(
    mcp_absolute_infil,
    "cell_type",
    "value",
    "group",
    my_y_positions = my_pvalues[["my_y_positions"]],
    my_x_min = my_pvalues[["my_x_min"]],
    my_x_max = my_pvalues[["my_x_max"]],
    my_annotation = my_pvalues[["my_annotation"]],
    my_palette,
    log_y = TRUE,
    cell_type = TRUE
  )

  # # transfor and plot per relative cell abundances
  # mcp_relative_infil = plot_tables[["mcp_relative_infil"]]
  #
  # # t-test
  # tests_relcell = mcp_relative_infil %>%
  #   dplyr::group_by(cell_type) %>%
  #   rstatix::t_test(value ~ group)
  # tests_relcell = tests_relcell %>% rstatix::add_xy_position(x = "cell_type")
  #
  # my_pvalues = get_pvalues_positions(mcp_relative_infil, tests_relcell)
  #
  #
  # #### CIBERSORT ####
  #
  # cb_all = plot_tables[["cb_all"]]
  #
  # ## retrieve sign annotation and sign positions
  # # t-test
  # tests_cb_all = cb_all %>%
  #   dplyr::group_by(cell_type) %>%
  #   rstatix::t_test(value ~ group)
  # tests_cb_all = tests_cb_all %>% rstatix::add_xy_position(x = "cell_type")
  #
  # my_pvalues = get_pvalues_positions(cb_all, tests_cb_all)
  #
  #
  # cb_sign = plot_tables[["cb_sign"]]
  #
  # ## retrieve sign annotation and sign positions
  # # t-test
  # tests_cb_sign = cb_sign %>%
  #   dplyr::group_by(cell_type) %>%
  #   rstatix::t_test(value ~ group)
  # tests_cb_sign = tests_cb_sign %>% rstatix::add_xy_position(x = "cell_type")
  #
  # my_pvalues = get_pvalues_positions(cb_sign, tests_cb_sign)


  #### quanTIseq ####

  mcp_relative_infil = plot_tables[["mcp_relative_infil"]]
  cb_all = plot_tables[["cb_all"]]
  quantiseq = plot_tables[["quantiseq"]]

  mcp_relative_infil$model = "MCP-counter"
  cb_all$model = "CIBERSORT"
  quantiseq$model = "quanTIseq"
  quantiseq$cell_type = stringr::str_replace_all(quantiseq$cell_type, "\\.", " ")

  all_relatives = mcp_relative_infil %>%
    BiocGenerics::rbind(cb_all) %>%
    BiocGenerics::rbind(quantiseq)

  # p_all_relatives = all_relatives %>%
  #   dplyr::group_by(cell_type) %>%
  #   rstatix::t_test(value ~ group, p.adjust.method = "BH")
  # p_all_relatives = p_all_relatives %>% rstatix::add_xy_position(x = "model")
  # p_all_relatives = p_all_relatives %>% rstatix::add_xy_position(x = "cell_type")

  p_all = ggplot(all_relatives, aes(x = cell_type, y = value, fill = group)) +
    geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
    stat_summary(fun = "median", colour = "red", geom = "crossbar", size = 0.15,
                 position = position_dodge(0.9), width = 0.5) +
    scale_fill_manual(values = my_palette) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    labs(x = "", y = "cell fractions", title = "Relative quantification") +
    facet_wrap(~ model, nrow = 3, scales = "free_x", strip.position = "right") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#d6eaf8"),
          strip.text = element_text(size = 10))


  # plot composite output figure
  p_infiltr = cowplot::plot_grid(
    cowplot::plot_grid(
      p_mcp_abs_all, p_mcp_abs_cells,
      nrow = 1,
      align = "h",
      rel_widths = c(0.3, 0.65)
    ),
    #pp_value,
    p_all,
    nrow = 2,
    align = "v",
    rel_heights = c(2, 5)
  )


  p_infiltr

  ##### TODO: add p-values for relative quantifications
  ##### TODO: add option to save figure png and pdf
  ##### TODO: return ggplot objects??


}
