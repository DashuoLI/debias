# ------------------------------
# LDA Utilities (complete)
# ------------------------------
library(MASS)

# Orthogonalize a vector `vec` against a list/matrix `basis`
orthogonalize <- function(vec, basis) {
  if (is.null(basis) || ncol(basis) == 0) return(vec / sqrt(sum(vec^2)))
  for (j in seq_len(ncol(basis))) {
    b <- basis[, j]
    vec <- vec - (sum(vec * b) / sum(b * b)) * b
  }
  vec / sqrt(sum(vec^2))
}

# Gram-Schmidt orthonormalization (columns of X)
gram_schmidt <- function(X, tol = 1e-12) {
  ncolX <- ncol(X)
  Q <- matrix(0, nrow = nrow(X), ncol = ncolX)
  k <- 0
  for (j in seq_len(ncolX)) {
    v <- X[, j]
    if (k > 0) {
      for (i in seq_len(k)) v <- v - sum(Q[, i] * X[, j]) * Q[, i]
    }
    norm_v <- sqrt(sum(v^2))
    if (norm_v > tol) {
      k <- k + 1
      Q[, k] <- v / norm_v
    }
  }
  if (k == 0) return(matrix(nrow = nrow(X), ncol = 0))
  Q[, seq_len(k), drop = FALSE]
}

# Complete basis: expand U to full orthonormal basis of R^N
complete_basis <- function(U) {
  N <- nrow(U)
  X <- cbind(U, diag(N))
  Q <- gram_schmidt(X)
  # Ensure exactly N columns (if possible)
  if (ncol(Q) < N) stop("Unable to build full basis (rank too low).")
  Q[, seq_len(N), drop = FALSE]
}

# Iterative LDA to extract orthogonal discriminant directions
perform_iterative_lda <- function(separating_pc, groups) {
  # separating_pc: samples x features (matrix)
  # groups: vector of group labels length = nrow(separating_pc)
  if (nrow(separating_pc) != length(groups)) stop("nrow(separating_pc) != length(groups)")
  n <- ncol(separating_pc)
  SCALING <- matrix(NA_real_, nrow = n, ncol = n)
  residuals <- separating_pc
  last_success <- 0

  for (idir in seq_len(n)) {
    lda_result <- tryCatch({
      MASS::lda(groups ~ ., data = as.data.frame(residuals))
    }, error = function(e) {
      message("lda() failed at iteration ", idir, ": ", e$message)
      NULL
    })

    if (is.null(lda_result)) break

    scaling <- as.numeric(lda_result$scaling)
    scaling <- scaling / sqrt(sum(scaling^2))
    if (idir > 1 && last_success >= 1) {
      # orthogonalize w.r.t previous scalings
      prev_basis <- if (last_success > 1) SCALING[, 1:last_success, drop = FALSE] else matrix(SCALING[, 1], ncol = 1)
      scaling <- orthogonalize(scaling, prev_basis)
    }

    # project residuals onto scaling and remove projection
    lda_x <- as.matrix(residuals) %*% scaling
    projection <- lda_x %*% t(scaling)
    residuals <- residuals - projection

    SCALING[, idir] <- scaling
    last_success <- idir
    message("finished LDA iteration ", idir)
  }

  if (last_success == 0) {
    stop("Failed LDA: no linear discriminant found.")
  }

  SCALING[, seq_len(last_success), drop = FALSE]
}

# adjust_rotation_with_lda: full wrapper similar to your original function
adjust_rotation_with_lda <- function(pca, annot, pcs_to_remove, fullxx_scaled, cancer_signal_annot,
                                    lda_p_cutoff = 0.05, lda_es_cutoff = 0.5, race_col = "Black",
                                    sample_id_col = "sample_id_alias", criterion = "p",
                                    eligibal_colvar_min_cutoff = 1e-08) {
  # If nothing to remove or single PC, return original rotation
  if (length(pcs_to_remove) == 0) {
    return(list(rotation = pca$rotation,
                pcs_to_remove = integer(0),
                pcs_to_keep = colnames(pca$rotation),
                stats = NULL))
  }
  if (length(pcs_to_remove) == 1) {
    return(list(rotation = pca$rotation,
                pcs_to_remove = pcs_to_remove,
                pcs_to_keep = setdiff(colnames(pca$rotation), paste0("PC", pcs_to_remove)),
                stats = NULL))
  }

  # prepare separating PCs (scores) for LDA
  pcs_idx <- pcs_to_remove
  separating_pc <- pca$x[annot[[sample_id_col]], pcs_idx, drop = FALSE]
  colvar_pc <- apply(separating_pc, 2, var)
  pcs_to_remove_low_var <- pcs_idx[colvar_pc <= eligibal_colvar_min_cutoff]
  pcs_for_lda <- pcs_idx[colvar_pc > eligibal_colvar_min_cutoff]

  if (length(pcs_for_lda) == 0) {
    # no pcs suitable for LDA
    return(list(rotation = pca$rotation,
                pcs_to_remove = pcs_to_remove_low_var,
                pcs_to_keep = setdiff(colnames(pca$rotation), paste0("PC", pcs_to_remove_low_var)),
                stats = NULL))
  }

  separating_pc_for_lda <- pca$x[annot[[sample_id_col]], paste0("PC", pcs_for_lda), drop = FALSE]
  # perform iterative LDA
  lda_res <- perform_iterative_lda(separating_pc_for_lda, annot[[race_col]])
  message("Finished perform_iterative_lda...")

  # rotation composition:
  rotation_part1 <- pca$rotation[, !(colnames(pca$rotation) %in% paste0("PC", pcs_to_remove)), drop = FALSE]
  rotation_part2 <- pca$rotation[, paste0("PC", pcs_for_lda), drop = FALSE] %*% lda_res
  rotation_part3 <- if (length(pcs_to_remove_low_var) > 0) pca$rotation[, paste0("PC", pcs_to_remove_low_var), drop = FALSE] else NULL

  # set column names
  colnames(rotation_part2) <- paste0("LDA", seq_len(ncol(rotation_part2)))
  if (!is.null(rotation_part3)) colnames(rotation_part3) <- paste0("LDAPC", seq_len(ncol(rotation_part3)))

  rotation <- cbind(rotation_part1, rotation_part2, rotation_part3)

  # compute new PCA-like scores after rotation for subsets and full data
  subsetxx_scaled_pca_scores_after_lda <- as.matrix(fullxx_scaled[annot[[sample_id_col]], , drop = FALSE]) %*% rotation
  fullxx_scaled_pca_scores_after_lda <- as.matrix(fullxx_scaled) %*% rotation

  # compute stats on new rotated scores (use compute_pc_stats from PCA utils)
  # compute_pc_stats expects pca-style inputs; we re-use that function by providing matrices directly
  # Note: compute_pc_stats defined in scripts/03_pca_utils.R should be available
  lda_stats <- compute_pc_stats(subsetxx_scaled_pca_scores_after_lda,
                                annot,
                                fullxx_scaled_pca_scores_after_lda,
                                cancer_signal_annot)

  # decide which LDA components to remove after tests (select_pcs_after_lda defined in pipeline)
  pcs_to_remove_after_lda <- select_pcs_after_lda(lda_stats, criterion = criterion,
                                                  padj_cutoff = lda_p_cutoff,
                                                  pblack_cutoff = 0, es_race_cutoff = lda_es_cutoff,
                                                  es_black_cutoff = 500.0)

  pcs_to_keep <- setdiff(lda_stats$pc_id, pcs_to_remove_after_lda)
  list(rotation = rotation,
       pcs_to_remove = pcs_to_remove_after_lda,
       pcs_to_keep = pcs_to_keep,
       stats = lda_stats)
}



select_pcs_after_lda <- function(pc_stats, criterion = "p", padj_cutoff = 0.05, pblack_cutoff = 0, es_race_cutoff = 0.5, es_black_cutoff = 500.0) {
  if (criterion == "p") {
    pc_stats$pc_id[pc_stats$padj < padj_cutoff & pc_stats$p_black_adj > pblack_cutoff & str_detect(pc_stats$pc_name, "LDA")]
  } else if (criterion == "es") {
    pc_stats$pc_id[abs(pc_stats$es_race) > es_race_cutoff & abs(pc_stats$es_black) < es_black_cutoff & str_detect(pc_stats$pc_name, "LDA")]
  }
  
}


