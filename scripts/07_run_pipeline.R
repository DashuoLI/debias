# ------------------------------
# Main pipeline wrapper
# ------------------------------

library(data.table)
library(MASS)
library(coin)
library(effectsize)

source("01_load_data.R")
source("02_stats_utils.R")
source("03_pca_utils.R")
source("04_lda_utils.R")
source("05_reconstruction_utils.R")
source("06_evaluation_utils.R")

run_pipeline <- function(marker_file, marker_type, annot_file, split,
                         data_file,
                         padj_cutoff = 0.05, pblack_cutoff = 0.01, lda_p_cutoff = 0.05,
                         select_pcs_criterion = "p", es_race_cutoff = 0.5, 
                         es_black_cutoff = 0.5, lda_es_cutoff = 0.5,
                         signal_test_set = "minor",
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
  output_file <- paste0("adjusted_", marker_type, 
                        ".select_pcs_criterion_", select_pcs_criterion, 
                        "__es_race_", es_race_cutoff,
                        "__es_black_", es_black_cutoff,
                        "__es_lda_", lda_es_cutoff,
                        "__signal_test_set_", signal_test_set,
                        ".csv.gz")
  fwrite(reconstructed_data_frame, output_file, row.names = F, col.names = T, sep = ",", quote = F, na = "NA")
  cat("adjusted data matrix written to:", output_file, "\n")
  message("Finished data output...")

  # Evaluate
  plot_file <- sub("\\.csv\\.gz$", ".stats.txt", output_file)
  eval_fn(noncancer_annot, mat_ori = fullxx_scaled,
          mat_adj = reconstructed, mat_min = fullxx_scaled - reconstructed,
          output_file = plot_file)
#  if (! is.null(lda_res$stats) ){
#    fwrite(lda_res$stats, paste(plot_file, "__lda_stats", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F, na = "NA")
#  }
#  fwrite(pc_stats, paste(plot_file, "__pc_stats", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F, na = "NA")
  cat("stats written in:", plot_file, "\n")
  message("Finished evaluation...")
}



