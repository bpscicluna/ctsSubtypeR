#' Classify new patients into molecular subtypes using packaged reference data
#'
#' Uses bundled reference gene expression data and subtype labels to classify
#' new samples into consensus transcriptomic subtypes (CTS). CTS group 4 is
#' treated as a health/null subtype.
#'
#' @param new_expr_data Numeric expression matrix or data frame of new samples
#'   with genes in rows and samples in columns. If a data frame, it may contain
#'   a first column named `gene` holding gene IDs.
#' @param gene_list Optional character vector of classifier genes. If `NULL`,
#'   all genes in the packaged reference dataset are used.
#' @param make_heatmap Logical; whether to draw a heatmap of the classified
#'   new samples. Default is `TRUE`.
#' @param ntrees Number of trees for random forest. Default is `500`.
#'
#' @return A list containing predicted subtypes, corrected expression matrices,
#'   silhouette metrics, and the trained random forest model.
#' @export
run_subtype_classifier <- function(new_expr_data,
                                   gene_list = NULL,
                                   make_heatmap = TRUE,
                                   ntrees = 500) {

  data("exp_core_g", package = "ctsSubtypeR", envir = environment())
  data("core_samples", package = "ctsSubtypeR", envir = environment())

  if (!exists("exp_core_g") || !exists("core_samples")) {
    stop("Packaged datasets `exp_core_g` and/or `core_samples` not found.")
  }

  if (!is.numeric(ntrees) || length(ntrees) != 1 || is.na(ntrees) || ntrees < 1) {
    stop("`ntrees` must be a positive integer.")
  }
  ntrees <- as.integer(ntrees)

  if (is.null(rownames(exp_core_g)) || anyDuplicated(rownames(exp_core_g))) {
    stop("`exp_core_g` must have unique gene IDs as row names.")
  }
  if (is.null(colnames(exp_core_g)) || anyDuplicated(colnames(exp_core_g))) {
    stop("`exp_core_g` must have unique sample IDs as column names.")
  }
  if (!is.data.frame(core_samples)) {
    stop("`core_samples` must be a data.frame.")
  }
  if (!"CTS" %in% colnames(core_samples)) {
    stop("`core_samples` must contain a column named `CTS`.")
  }
  if (is.null(rownames(core_samples)) || anyDuplicated(rownames(core_samples))) {
    stop("`core_samples` must have unique row names matching `colnames(exp_core_g)`.")
  }

  missing_core_meta <- setdiff(colnames(exp_core_g), rownames(core_samples))
  if (length(missing_core_meta) > 0) {
    stop(
      "These core samples are missing in `core_samples`: ",
      paste(missing_core_meta, collapse = ", ")
    )
  }
  core_samples <- core_samples[colnames(exp_core_g), , drop = FALSE]

  cts <- as.factor(core_samples$CTS)
  if (any(is.na(cts))) {
    stop("`core_samples$CTS` contains missing values.")
  }

  valid_cts <- c("1", "2", "3", "4")
  if (!all(as.character(cts) %in% valid_cts)) {
    stop("`core_samples$CTS` must only contain: 1, 2, 3, or 4.")
  }

  cts <- droplevels(factor(as.character(cts)))
  train_levels <- levels(cts)

  if (is.data.frame(new_expr_data)) {
    if ("gene" %in% colnames(new_expr_data)) {
      rownames(new_expr_data) <- new_expr_data$gene
      new_expr_data <- new_expr_data[, setdiff(colnames(new_expr_data), "gene"), drop = FALSE]
    } else if (!is.numeric(new_expr_data[[1]])) {
      rownames(new_expr_data) <- new_expr_data[[1]]
      new_expr_data <- new_expr_data[, -1, drop = FALSE]
    }
    new_expr_data <- as.matrix(new_expr_data)
  }

  if (!is.matrix(new_expr_data)) {
    stop("`new_expr_data` must be a matrix or data.frame.")
  }
  if (is.null(rownames(new_expr_data)) || anyDuplicated(rownames(new_expr_data))) {
    stop("`new_expr_data` must have unique gene IDs as row names, or a gene column.")
  }
  if (is.null(colnames(new_expr_data)) || anyDuplicated(colnames(new_expr_data))) {
    stop("`new_expr_data` must have unique sample IDs as column names.")
  }

  storage.mode(exp_core_g) <- "numeric"
  storage.mode(new_expr_data) <- "numeric"

  if (is.null(gene_list)) {
    gene_list <- rownames(exp_core_g)
  } else {
    gene_list <- unique(as.character(gene_list))
  }

  genes_use <- intersect(gene_list, rownames(exp_core_g))
  genes_use <- intersect(genes_use, rownames(new_expr_data))

  if (length(genes_use) < 10) {
    stop("Too few overlapping genes between reference and new data: ", length(genes_use), ".")
  }

  exp_core_sub <- exp_core_g[genes_use, , drop = FALSE]
  new_data_sub <- new_expr_data[genes_use, , drop = FALSE]

  exp_all_m <- cbind(exp_core_sub, new_data_sub)
  batch <- c(rep("core", ncol(exp_core_sub)), rep("new", ncol(new_data_sub)))

  exp_all_combat <- sva::ComBat(
    dat = exp_all_m,
    batch = batch,
    par.prior = TRUE,
    prior.plots = FALSE
  )

  exp_core_combat <- exp_all_combat[, colnames(exp_core_sub), drop = FALSE]
  exp_new_combat <- exp_all_combat[, colnames(new_data_sub), drop = FALSE]

  train_df <- as.data.frame(t(exp_core_combat))
  train_df$CTS <- cts

  rf <- randomForest::randomForest(
    CTS ~ .,
    data = train_df,
    ntree = ntrees,
    proximity = TRUE
  )

  pred_new <- stats::predict(rf, newdata = as.data.frame(t(exp_new_combat)))
  pred_new <- factor(as.character(pred_new), levels = train_levels)

  pred_df <- data.frame(
    sample_id = colnames(exp_new_combat),
    CTS = pred_new
  )

  if (isTRUE(make_heatmap)) {
    ordered_pred_df <- pred_df[order(pred_df$CTS), , drop = FALSE]
    rownames(ordered_pred_df) <- ordered_pred_df$sample_id
    exp_new_ordered <- exp_new_combat[, ordered_pred_df$sample_id, drop = FALSE]

    ann_colors <- list(
      CTS = c("1" = "royalblue", "2" = "#B2DF8A", "3" = "orange", "4" = "grey70")
    )

    pheatmap::pheatmap(
      exp_new_ordered,
      color = grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(11, "PuOr"))
      )(50),
      clustering_distance_rows = "correlation",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = FALSE,
      annotation_col = ordered_pred_df[, "CTS", drop = FALSE],
      annotation_colors = ann_colors,
      scale = "row",
      treeheight_row = 0,
      border_color = NA
    )
  }

  sil <- NULL
  keep <- !is.na(pred_new)
  if (sum(keep) > 2 && length(unique(pred_new[keep])) > 1) {
    rf_new <- randomForest::randomForest(
      x = as.data.frame(t(exp_new_combat[, keep, drop = FALSE])),
      y = droplevels(pred_new[keep]),
      ntree = ntrees,
      proximity = TRUE
    )

    sil <- cluster::silhouette(
      as.integer(droplevels(pred_new[keep])),
      stats::as.dist(1 - rf_new$proximity)
    )
  }

  list(
    predictions = pred_df,
    expression_corrected = list(core = exp_core_combat, new_data = exp_new_combat),
    silhouette = sil,
    rf_model = rf,
    genes_used = genes_use,
    gene_overlap = length(genes_use)
  )
}
