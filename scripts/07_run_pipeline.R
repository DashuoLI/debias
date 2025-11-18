# ------------------------------
# Main pipeline wrapper
# ------------------------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(effectsize))

source("01_load_data.R")
source("02_stats_utils.R")
source("03_pca_utils.R")
source("04_lda_utils.R")
source("05_reconstruction_utils.R")
source("06_evaluation_utils.R")

run_pipeline <- function(marker_file, annot_file,
                         data_file, output_dir,
                         padj_cutoff = 0.05, pblack_cutoff = 0.05, lda_p_cutoff = 0.05,
                         select_pcs_criterion = "p", es_race_cutoff = 0.1, 
                         es_black_cutoff = 0.5, lda_es_cutoff = 0.1,
                         signal_test_set = "full",
                         eval_fn = evaluate_features, 
                         sample_id_col = "sample_id_alias", 
                         race_col = "Black") {
  
  noncancer_annot <- load_annotation(annot_file, filter_label = "both")
  demographic_annot <- load_annotation(annot_file, filter_label = "demographic")
  full_annot <- load_annotation(annot_file, filter_label = "none")  
  dt <- load_data(data_file)
  markers <- load_markers(marker_file, dt)
  mat <- as.matrix(dt[markers, noncancer_annot[[sample_id_col]]])
  message("Finished data loading...")

  # PCA
  pca_res <- perform_pca(mat)
  pca <- pca_res$pca
  colvar <- pca_res$colvar
  message("Finished PCA...")
 
  # Full matrix scaling
  fullmat <- as.matrix(dt[markers, -1])
  fullxx <- t(fullmat)
  fullxx_zerovar <- fullxx[, colvar <= eligibal_colvar_min_cutoff, drop = FALSE]
  fullxx <- fullxx[, colvar > eligibal_colvar_min_cutoff, drop = FALSE]
  fullxx_scaled <- scale(fullxx, center = pca$center, scale = pca$scale)
  fullxx_scaled_pca_scores <- as.matrix(fullxx_scaled) %*% pca$rotation
  message("Finished full matrix scaling...")

  # Biased PC selection
  if (signal_test_set == "minor") {
    cancer_signal_annot <- demographic_annot[demographic_annot[[race_col]] == "Black", ]
  } else if (signal_test_set == "full") {
    cancer_signal_annot <- demographic_annot
  }
  pc_stats <- compute_pc_stats(pca$x, noncancer_annot, fullxx_scaled_pca_scores, cancer_signal_annot)
  pcs_to_remove <- select_pcs(pc_stats, criterion = select_pcs_criterion, padj_cutoff = padj_cutoff, pblack_cutoff = pblack_cutoff, es_race_cutoff = es_race_cutoff, es_black_cutoff = es_black_cutoff)
  message("Finished biased PC selection...")
  cat("pcs_to_remove_before_lda:", pcs_to_remove, "\n")

  # LDA adjustment
  lda_res <- adjust_rotation_with_lda(pca, noncancer_annot, pcs_to_remove, fullxx_scaled, cancer_signal_annot,  lda_p_cutoff, lda_es_cutoff, criterion = select_pcs_criterion)
  message("Finished LDA adjustment...")
  lda_stats <- lda_res$stats
  cat("pcs_to_remove_after_lda:", lda_stats[lda_res$pcs_to_remove, "pc_name"], "\n")

  # Reconstruction
  reconstructed <- reconstruct_matrix(fullxx_scaled, lda_res$rotation,
                                      lda_res$pcs_to_keep)
  reconstructed_nonneg <- make_nonnegative(reconstructed, full_annot, fullxx_zerovar)
  reconstructed_data_frame <- make_output_data_frame(reconstructed_nonneg)
  message("Finished reconstruction...")

  # Save
  # Check if output_dir exists, and create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)  # recursive = TRUE creates nested directories
  }
  output_file <- paste0(output_dir, "/", "adjusted_data.csv.gz")
  fwrite(reconstructed_data_frame, output_file, row.names = F, col.names = T, sep = ",", quote = F, na = "NA")
  cat("adjusted data matrix written to:", output_file, "\n")
  message("Finished data output...")

}



