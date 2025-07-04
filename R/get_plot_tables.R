### Prep estimates tables for plotting -----------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for plotting infiltR summary figures
#' @note 2024-11-13
#' @title get_plot_tables
#' @details
#'
#' This function prepr estimates tables for plotting
#'
#' @param infiltr_out infiltR output object
#' @param metadata A metadata matrix with samples as rows
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
#' @param sample_levels Order of samples groups for plotting.
#'   Default: "default", groups are plotted in alphabetical order.
#' @export
get_plot_tables = function(
    infiltr_out,
    metadata,
    sample_groups,
    sample_levels
){

  # declare output obj
  plot_tables = list()


  #### MCP-counter ####

  mcp_absolute_infil = infiltr_out[["mcp_counter"]]

  # transfor and plot per cell abundances
  mcp_absolute_infil$Sample = rownames(mcp_absolute_infil)
  mcp_absolute_infil$group = metadata[, sample_groups]
  mcp_absolute_infil = mcp_absolute_infil %>%
    dplyr::select(-c("all_infiltrating"))
  mcp_absolute_infil = reshape2::melt(mcp_absolute_infil)
  colnames(mcp_absolute_infil) = c("group", "Sample", "cell_type", "value")

  mcp_absolute_infil$cell_type = factor(
    mcp_absolute_infil$cell_type,
    levels = c(
      "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage",
      "NK cells", "Monocytic lineage", "Myeloid dendritic cells",
      "Neutrophils", "Endothelial cells", "Fibroblasts"
    )
  )

  # transfor and plot per relative cell abundances
  mcp_relative_infil = infiltr_out[["mcp_counter_norm"]]
  mcp_relative_infil$Sample = rownames(mcp_relative_infil)
  mcp_relative_infil$group = metadata[, sample_groups]
  mcp_relative_infil = reshape2::melt(mcp_relative_infil)
  mcp_relative_infil$value = mcp_relative_infil$value %>%
    as.numeric() %>%
    round(digits = 4)
  colnames(mcp_relative_infil) = c("group", "Sample", "cell_type", "value")

  mcp_relative_infil$cell_type = factor(
    mcp_relative_infil$cell_type,
    levels = c(
      "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage",
      "NK cells", "Monocytic lineage", "Myeloid dendritic cells",
      "Neutrophils", "Endothelial cells", "Fibroblasts"
    )
  )



  #### CIBERSORT ####

  cb_all = infiltr_out[["cb_all"]]

  # transfor and plot per cell abundances
  cb_all$Sample = rownames(cb_all)
  cb_all$group = metadata[, sample_groups]
  cb_all = cb_all %>%
    dplyr::select(-c("P-value", "Correlation", "RMSE"))
  cb_all = reshape2::melt(cb_all)
  colnames(cb_all) = c("Sample", "group", "cell_type", "value")

  cb_all$cell_type = factor(
    cb_all$cell_type,
    levels = c(
      "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
      "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
      "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta",
      "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0",
      "Macrophages M1", "Macrophages M2", "Dendritic cells resting",
      "Dendritic cells activated", "Mast cells resting", "Mast cells activated",
      "Eosinophils", "Neutrophils"
    )
  )

  # get significant CIBERSORT samples if not empty
  cb_sign = infiltr_out[["cb_sign"]]

  if(nrow(cb_sign) >= 1){

    # transform and plot per cell abundances
    cb_sign$Sample = rownames(cb_sign)
    cb_sign$group = metadata %>%
      dplyr::filter(sample %in% rownames(cb_sign)) %>%
      dplyr::select(dplyr::all_of(sample_groups)) %>%
      BiocGenerics::unlist()
    cb_sign = cb_sign %>%
      dplyr::select(-c("P-value", "Correlation", "RMSE"))
    cb_sign = reshape2::melt(cb_sign)
    colnames(cb_sign) = c("Sample", "group", "cell_type", "value")

    cb_sign$cell_type = factor(
      cb_sign$cell_type,
      levels = c(
        "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
        "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
        "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta",
        "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0",
        "Macrophages M1", "Macrophages M2", "Dendritic cells resting",
        "Dendritic cells activated", "Mast cells resting", "Mast cells activated",
        "Eosinophils", "Neutrophils"
      )
    )

  }



  #### quanTIseq ####

  quantiseq = infiltr_out[["quantiseq"]]

  # transfor and plot per cell abundances
  quantiseq$Sample = rownames(quantiseq)
  quantiseq$group = metadata[, sample_groups]
  quantiseq = reshape2::melt(quantiseq)
  colnames(quantiseq) = c("Sample", "group", "cell_type", "value")

  quantiseq$cell_type = factor(
    quantiseq$cell_type,
    levels = c(
      "B.cells", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils",
      "NK.cells", "T.cells.CD4", "T.cells.CD8", "Tregs", "Dendritic.cells", "Other"
    )
  )



  # order samples if requested
  mcp_absolute_infil$group = factor(mcp_absolute_infil$group, levels = sample_levels)
  mcp_relative_infil$group = factor(mcp_relative_infil$group, levels = sample_levels)
  cb_all$group = factor(cb_all$group, levels = sample_levels)
  cb_sign$group = factor(cb_sign$group, levels = sample_levels)
  quantiseq$group = factor(quantiseq$group, levels = sample_levels)

  # add tables
  plot_tables[["mcp_absolute_infil"]] = mcp_absolute_infil
  plot_tables[["mcp_relative_infil"]] = mcp_relative_infil
  plot_tables[["cb_all"]] = cb_all
  plot_tables[["cb_sign"]] = cb_sign
  plot_tables[["quantiseq"]] = quantiseq

  return(plot_tables)

}
