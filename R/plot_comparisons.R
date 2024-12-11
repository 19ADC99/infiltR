### Plot estimates comparisons -------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note plot estimates comparisons
#' @note 2024-12-09
#'
#' @param infiltr_out infiltr_out output object
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
plot_comparisons = function(
    infiltr_out,
    save_plots = TRUE,
    outdir = "default"
){

  #### MCP-counter VS CB all ####
  mcp_counter_norm = infiltr_out[["mcp_counter_norm"]]
  cb_all = infiltr_out[["cb_all"]]

  mcp_counter_norm$Sample = rownames(mcp_counter_norm)
  cb_all$Sample = rownames(cb_all)

  cells_to_keep = c(
    "Sample", "T.cells", "T.cells.CD8", "B.cells", "NK.cells",
    "Monocytes", "Dendritic.cells", "Neutrophils"
  )

  ## merge CIBERSORT cell types to match MCP-counter
  cb_all$B.cells = cb_all[, "B cells naive"] + cb_all[, "B cells memory"]
  cb_all$NK.cells = cb_all[, "NK cells resting"] + cb_all[, "NK cells activated"]
  cb_all$Dendritic.cells = cb_all[, "Dendritic cells resting"] + cb_all[, "Dendritic cells activated"]
  cb_all$T.cells = cb_all[, "T cells CD4 naive"] +
    cb_all[, "T cells CD4 memory resting"] +
    cb_all[, "T cells CD4 memory activated"] +
    cb_all[, "T cells follicular helper"] +
    cb_all[, "T cells regulatory (Tregs)"] +
    cb_all[, "T cells gamma delta"]
  cb_all$Monocytes =  cb_all[, "Monocytes"] +
    cb_all[, "Macrophages M0"] +
    cb_all[, "Macrophages M1"] +
    cb_all[, "Macrophages M2"]

  cb_all = cb_all %>%
    dplyr::rename("T.cells.CD8" = "T cells CD8") %>%
    dplyr::select(all_of(cells_to_keep)) %>%
    data.table::melt()
  colnames(cb_all) = c("Sample", "cell_type", "est_cb")

  mcp_counter_norm = mcp_counter_norm %>%
    dplyr::rename(
      "T.cells.CD8" = "CD8 T cells",
      "T.cells" = "T cells",
      "B.cells" = "B lineage",
      "NK.cells" = "NK cells",
      "Monocytes" = "Monocytic lineage",
      "Dendritic.cells" = "Myeloid dendritic cells"
    ) %>%
    dplyr::select(all_of(cells_to_keep)) %>%
    data.table::melt()
  colnames(mcp_counter_norm) = c("Sample", "cell_type", "est_mcp")

  compare_df = mcp_counter_norm
  compare_df$est_cb = cb_all$est_cb

  p_mcp_vs_cb = ggplot(compare_df, aes(x = est_mcp, y = est_cb)) +
    geom_point(size = 0.25, color = "grey55") +
    geom_smooth(method = lm, se = FALSE, color = "firebrick") +
    ggpmisc::stat_correlation() +
    facet_wrap(~ cell_type, ncol = 5, scales = "free") +
    labs(x = "MCP-counter", y = "CIBERSORT", fill = "") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#efd2d2"),
          strip.text = element_text(size = 10))



  #### quanTIseq VS CB all ####
  quantiseq = infiltr_out[["quantiseq"]]
  cb_all = infiltr_out[["cb_all"]]

  cb_all$Sample = rownames(cb_all)

  cells_to_keep = c(
    "Sample", "B.cells", "Macrophages.M1", "Macrophages.M2", "Monocytes",
    "Neutrophils", "NK.cells", "T.cells.CD4", "T.cells.CD8", "Tregs", "Dendritic.cells"
  )

  ## merge CIBERSORT cell types to match quanTIseq
  cb_all$B.cells = cb_all[, "B cells naive"] + cb_all[, "B cells memory"]
  cb_all$NK.cells = cb_all[, "NK cells resting"] + cb_all[, "NK cells activated"]
  cb_all$Dendritic.cells = cb_all[, "Dendritic cells resting"] + cb_all[, "Dendritic cells activated"]
  cb_all$T.cells.CD4 = cb_all[, "T cells CD4 naive"] +
    cb_all[, "T cells CD4 memory resting"] +
    cb_all[, "T cells CD4 memory activated"] +
    cb_all[, "T cells follicular helper"]

  cb_all = cb_all %>%
    dplyr::rename(
      "Macrophages.M1" = "Macrophages M1",
      "Macrophages.M2" = "Macrophages M2",
      "T.cells.CD8" = "T cells CD8",
      "Tregs" = "T cells regulatory (Tregs)"
    ) %>%
    dplyr::select(all_of(cells_to_keep)) %>%
    data.table::melt()
  colnames(cb_all) = c("Sample", "cell_type", "est_cb")

  quantiseq = quantiseq %>%
    dplyr::select(all_of(cells_to_keep))%>%
    data.table::melt()
  colnames(quantiseq) = c("Sample", "cell_type", "est_quantiseq")

  compare_df = cb_all
  compare_df$est_quantiseq = quantiseq$est_quantiseq

  p_cb_vs_quantiseq = ggplot(compare_df, aes(x = est_cb, y = est_quantiseq)) +
    geom_point(size = 0.25, color = "grey55") +
    geom_smooth(method = lm, se = FALSE, color = "steelblue") +
    ggpmisc::stat_correlation() +
    facet_wrap(~ cell_type, ncol = 5, scales = "free") +
    labs(x = "CIBERSORT", y = "quanTIseq", fill = "") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#d6eaf8"),
          strip.text = element_text(size = 10))



  #### quanTIseq VS MCP-counter ####
  mcp_counter_norm = infiltr_out[["mcp_counter_norm"]]
  quantiseq = infiltr_out[["quantiseq"]]

  mcp_counter_norm$Sample = rownames(mcp_counter_norm)

  cells_to_keep = c(
    "Sample", "B.cells", "Monocytes", "Neutrophils", "NK.cells",
    "T.cells.CD8", "T.cells", "Dendritic.cells"
  )

  mcp_counter_norm = mcp_counter_norm %>%
    dplyr::rename(
      "T.cells.CD8" = "CD8 T cells",
      "T.cells" = "T cells",
      "B.cells" = "B lineage",
      "NK.cells" = "NK cells",
      "Monocytes" = "Monocytic lineage",
      "Dendritic.cells" = "Myeloid dendritic cells"
    ) %>%
    dplyr::select(all_of(cells_to_keep)) %>%
    data.table::melt()
  colnames(mcp_counter_norm) = c("Sample", "cell_type", "est_mcp")

  quantiseq$T.cells = quantiseq[, "T.cells.CD4"] + quantiseq[, "Tregs"]
  quantiseq$Monocytes = quantiseq[, "Monocytes"] +
    quantiseq[, "Macrophages.M1"] +
    quantiseq[, "Macrophages.M2"]
  quantiseq = quantiseq %>%
    dplyr::select(all_of(cells_to_keep))%>%
    data.table::melt()
  colnames(quantiseq) = c("Sample", "cell_type", "est_quantiseq")

  compare_df = mcp_counter_norm
  compare_df$est_quantiseq = quantiseq$est_quantiseq

  p_mcp_vs_quantiseq = ggplot(compare_df, aes(x = est_mcp, y = est_quantiseq)) +
    geom_point(size = 0.25, color = "grey55") +
    geom_smooth(method = lm, se = FALSE, color = "orange") +
    ggpmisc::stat_correlation() +
    facet_wrap(~ cell_type, ncol = 5, scales = "free") +
    labs(x = "MCP-counter", y = "quanTIseq", fill = "") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#ffd7b5"),
          strip.text = element_text(size = 10))


  # plot composite output figure
  p_comparisons = cowplot::plot_grid(
    p_mcp_vs_cb,
    NULL,
    p_cb_vs_quantiseq,
    NULL,
    p_mcp_vs_quantiseq,
    nrow = 5,
    rel_heights = c(1, 0.1, 1, 0.1, 1)
  )

  p_comparisons

}
