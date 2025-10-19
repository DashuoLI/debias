# ------------------------------
# PCA Utilities
# ------------------------------

eligibal_colvar_min_cutoff <- 1.0e-08

perform_pca <- function(mat) {
  x <- t(mat)
  colvar <- apply(x, 2, var)
  x <- x[, colvar > eligibal_colvar_min_cutoff]
  pca <- prcomp(x, scale. = TRUE)
  list(pca = pca, colvar = colvar)
}

compute_pc_stats <- function(race_pca_scores, race_annot, cancer_pca_scores, cancer_annot, race_col = "Group", signal_col = "Label", sample_id_col = "sample_id_alias") {
  n_pc <- ncol(race_pca_scores)
  stats <- data.frame(pc_id = 1:n_pc, pc_name = colnames(race_pca_scores), p_race = NA_real_, p_black = NA_real_, es_race = NA_real_, es_black = NA_real_)
  
  for (i in 1:n_pc) {
    race_pc_values <- race_pca_scores[, i]
    race_annot_groups <- race_annot[[race_col]]
    stats$p_race[i] <- safe_fisher_pitman(race_pc_values, race_annot[[race_col]])
    stats$es_race[i] <- safe_cohens_d(race_pc_values[race_annot_groups == "Black-noncancer"], race_pc_values[race_annot_groups == "White-noncancer"], pooled_sd = TRUE)   
 
    cancer_pc_values <- cancer_pca_scores[cancer_annot[[sample_id_col]], i]
    cancer_annot_groups <- cancer_annot[[signal_col]]
    stats$p_black[i] <- safe_fisher_pitman(cancer_pc_values, cancer_annot_groups)
    stats$es_black[i] <- safe_cohens_d(cancer_pc_values[cancer_annot_groups == "cancer"], cancer_pc_values[cancer_annot_groups == "noncancer"], pooled_sd = TRUE)  
  }
  
  stats$padj <- p.adjust(stats$p_race, method = "BH") ## adjust for p_race and p_cancer
  stats$p_black_adj <- p.adjust(stats$p_black, method = "BH")
  stats
}

select_pcs <- function(pc_stats, criterion = "p", padj_cutoff = 0.05, pblack_cutoff = 0.05, es_race_cutoff = 0.5, es_black_cutoff = 0.5) {
  if (criterion == "p") {
    pc_stats$pc_id[pc_stats$padj < padj_cutoff & pc_stats$p_black_adj > pblack_cutoff]
  } else if (criterion == "es") {
    pc_stats$pc_id[abs(pc_stats$es_race) > es_race_cutoff & abs(pc_stats$es_black) < es_black_cutoff]
  }
}


