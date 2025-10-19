# ------------------------------
# Reconstruction Utilities
# ------------------------------

# reconstruct_matrix: projects full scaled data to rotated basis, zeroes out retained PCs, reconstructs
reconstruct_matrix <- function(fullxx_scaled, rotation, pcs_to_keep) {
  # fullxx_scaled: samples x features (centered & scaled)
  # rotation: features x components
  # pcs_to_keep: vector of pc indexes/names to set to zero (these are indices in pc space)
  scores <- as.matrix(fullxx_scaled) %*% rotation
  # pcs_to_keep can be numeric indices of columns in scores or names "LDA1" etc.
  if (is.character(pcs_to_keep)) {
    keep_idx <- which(colnames(scores) %in% pcs_to_keep)
  } else {
    keep_idx <- pcs_to_keep
  }
  if (length(keep_idx) > 0) {
    scores[, keep_idx] <- 0
  }
  reconstructed <- fullxx_scaled - scores %*% t(rotation)
  reconstructed
}

# Make result non-negative relative to sample minima (preserve shape)
make_nonnegative <- function(mat, annot, fullxx_zerovar, sample_id_col = "sample_id_alias") {
  # mat: reconstructed (samples x features)
  # annot: full annotation (contains sample IDs)
  # fullxx_zerovar: columns that were removed earlier due to zero var - appended as zeros
  # Guarantee non-negativity per sample across provided sample list
  if (!is.null(fullxx_zerovar) && ncol(fullxx_zerovar) > 0) {
    mat_full <- cbind(mat, fullxx_zerovar)
  } else {
    mat_full <- mat
  }
  # Shift each row so minimum is zero
  row_mins <- apply(mat_full, 1, min, na.rm = TRUE)
  t(t(mat_full) - row_mins)
}

make_output_data_frame <- function(reconstructed) {
  dt <- as.data.frame(reconstructed)
  dt$key <- rownames(reconstructed)
  # put key as first column
  dt <- dt[, c("key", setdiff(names(dt), "key"))]
  dt
}


