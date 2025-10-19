# ------------------------------
# Evaluation Utilities (full)
# ------------------------------

library(infotheo)

evaluate_features <- function(annot, mat_ori, mat_adj, mat_min,
                              group_col = "Black", group_num_col = "Black2",
                              output_file = NULL, sample_id_col = "sample_id_alias") {
  feature_names <- colnames(mat_adj)
  subset_mat_ori <- mat_ori[annot[[sample_id_col]],]
  subset_mat_adj <- mat_adj[annot[[sample_id_col]],]
  subset_mat_min <- mat_min[annot[[sample_id_col]],]

  stats <- data.frame(
    key = feature_names,
    tp_adj = NA_real_, tp_ori = NA_real_, tp_min = NA_real_,
    wp_adj = NA_real_, wp_ori = NA_real_, wp_min = NA_real_,
    fpp_adj = NA_real_, fpp_ori = NA_real_, fpp_min = NA_real_,
    es_adj = NA_real_, es_ori = NA_real_, es_min = NA_real_,
    mutinfo_adj = NA_real_, mutinfo_ori = NA_real_, mutinfo_min = NA_real_,
    shapiro_adj_blk = NA_real_, shapiro_adj_wht = NA_real_, 
    shapiro_ori_blk = NA_real_, shapiro_ori_wht = NA_real_
  )
  
  for (i in seq_along(feature_names)) {
    g_black <- annot[[group_col]] == "Black"
    g_white <- annot[[group_col]] == "White"
    
    stats$tp_adj[i] <- safe_ttest(subset_mat_adj[g_black, i], subset_mat_adj[g_white, i])
    stats$tp_ori[i] <- safe_ttest(subset_mat_ori[g_black, i], subset_mat_ori[g_white, i])
    stats$tp_min[i] <- safe_ttest(subset_mat_min[g_black, i], subset_mat_min[g_white, i])
    
    stats$wp_adj[i] <- safe_wilcox(subset_mat_adj[g_black, i], subset_mat_adj[g_white, i])
    stats$wp_ori[i] <- safe_wilcox(subset_mat_ori[g_black, i], subset_mat_ori[g_white, i])
    stats$wp_min[i] <- safe_wilcox(subset_mat_min[g_black, i], subset_mat_min[g_white, i])

    stats$fpp_adj[i] <- safe_fisher_pitman(subset_mat_adj[, i], annot[[group_col]])
    stats$fpp_ori[i] <- safe_fisher_pitman(subset_mat_ori[, i], annot[[group_col]])
    stats$fpp_min[i] <- safe_fisher_pitman(subset_mat_min[, i], annot[[group_col]])

    stats$shapiro_adj_blk[i] <- safe_shapiro(subset_mat_adj[g_black, i])
    stats$shapiro_adj_wht[i] <- safe_shapiro(subset_mat_adj[g_white, i])
    stats$shapiro_ori_blk[i] <- safe_shapiro(subset_mat_ori[g_black, i])
    stats$shapiro_ori_wht[i] <- safe_shapiro(subset_mat_ori[g_white, i]) 
    
    stats$mutinfo_adj[i] <- mutinformation(discretize(subset_mat_adj[, i]), annot[[group_col]], method = "emp")
    stats$mutinfo_ori[i] <- mutinformation(discretize(subset_mat_ori[, i]), annot[[group_col]], method = "emp")
    stats$mutinfo_min[i] <- mutinformation(discretize(subset_mat_min[, i]), annot[[group_col]], method = "emp")

    stats$es_adj[i] <- safe_cohens_d(subset_mat_adj[g_black, i], subset_mat_adj[g_white, i], pooled_sd = TRUE)
    stats$es_ori[i] <- safe_cohens_d(subset_mat_ori[g_black, i], subset_mat_ori[g_white, i], pooled_sd = TRUE)
    stats$es_min[i] <- safe_cohens_d(subset_mat_min[g_black, i], subset_mat_min[g_white, i], pooled_sd = TRUE)

  }
  
  if (!is.null(output_file)) {
    fwrite(stats, output_file, sep = "\t", quote = FALSE, na = "NA")
    cat("stats written in:", output_file, "\n")
  }
  stats
}



