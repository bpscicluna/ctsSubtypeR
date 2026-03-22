#' Normalize reference, new healthy, and new case samples with ComBat
#'
#' Performs joint batch correction across packaged reference samples,
#' user-supplied healthy controls, and user-supplied case samples while
#' preserving healthy-versus-sepsis biology.
#'
#' @param exp_core_g Reference expression matrix; genes in rows, samples in columns.
#' @param core_samples Reference sample annotation data frame containing `CTS`.
#' @param new_case_data Matrix or data frame of new case samples.
#' @param new_healthy_data Matrix or data frame of new healthy control samples.
#' @param gene_list Optional character vector of classifier genes.
#'
#' @return A list containing corrected matrices, genes used, batch labels, and biological group labels.
#' @keywords internal
normalize_with_healthy_combat <- function(exp_core_g,
                                          core_samples,
                                          new_case_data,
                                          new_healthy_data,
                                          gene_list = NULL) {

  if (is.null(rownames(exp_core_g))) {
    stop("`exp_core_g` must have row names with gene IDs.")
  }
  if (is.null(colnames(exp_core_g))) {
    stop("`exp_core_g` must have sample IDs in columns.")
  }
  if (!"CTS" %in% colnames(core_samples)) {
    stop("`core_samples` must contain a `CTS` column.")
  }
  if (is.null(rownames(core_samples))) {
    stop("`core_samples` must have row names matching `colnames(exp_core_g)`.")
  }

  core_samples <- core_samples[colnames(exp_core_g), , drop = FALSE]

  prep_expr <- function(x, object_name) {
    if (is.data.frame(x)) {
      if ("gene" %in% colnames(x)) {
        rownames(x) <- x$gene
        x <- x[, setdiff(colnames(x), "gene"), drop = FALSE]
      } else if (!is.numeric(x[[1]])) {
        rownames(x) <- x[[1]]
        x <- x[, -1, drop = FALSE]
      }
      x <- as.matrix(x)
    }

    if (!is.matrix(x)) {
      stop(sprintf("`%s` must be a matrix or data.frame.", object_name))
    }
    if (is.null(rownames(x))) {
      stop(sprintf("`%s` must have row names with gene IDs or a gene column.", object_name))
    }
    if (is.null(colnames(x))) {
      stop(sprintf("`%s` must have sample IDs in columns.", object_name))
    }

    storage.mode(x) <- "numeric"
    x
  }

  new_case_data <- prep_expr(new_case_data, "new_case_data")
  new_healthy_data <- prep_expr(new_healthy_data, "new_healthy_data")
  storage.mode(exp_core_g) <- "numeric"

  if (is.null(gene_list)) {
    gene_list <- rownames(exp_core_g)
  } else {
    gene_list <- unique(as.character(gene_list))
  }

  genes_use <- Reduce(
    intersect,
    list(
      gene_list,
      rownames(exp_core_g),
      rownames(new_case_data),
      rownames(new_healthy_data)
    )
  )

  if (length(genes_use) < 10) {
    stop("Too few overlapping genes across reference, new cases, and new healthy controls.")
  }

  exp_core_sub <- exp_core_g[genes_use, , drop = FALSE]
  case_sub <- new_case_data[genes_use, , drop = FALSE]
  healthy_sub <- new_healthy_data[genes_use, , drop = FALSE]

  exp_all <- cbind(exp_core_sub, healthy_sub, case_sub)

  batch <- c(
    rep("reference", ncol(exp_core_sub)),
    rep("new_healthy", ncol(healthy_sub)),
    rep("new_case", ncol(case_sub))
  )

  ref_group <- ifelse(as.character(core_samples$CTS) == "4", "healthy", "sepsis")

  bio_group <- c(
    ref_group,
    rep("healthy", ncol(healthy_sub)),
    rep("sepsis", ncol(case_sub))
  )

  batch <- factor(batch)
  bio_group <- factor(bio_group, levels = c("healthy", "sepsis"))

  mod <- stats::model.matrix(~ bio_group)

  exp_all_combat <- sva::ComBat(
    dat = exp_all,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )

  list(
    corrected = list(
      reference = exp_all_combat[, colnames(exp_core_sub), drop = FALSE],
      new_healthy = exp_all_combat[, colnames(healthy_sub), drop = FALSE],
      new_case = exp_all_combat[, colnames(case_sub), drop = FALSE]
    ),
    genes_used = genes_use,
    batch = batch,
    bio_group = bio_group
  )
}
