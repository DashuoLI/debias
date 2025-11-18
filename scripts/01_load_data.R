# ------------------------------
# Data Loading Utilities
# ------------------------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

load_annotation <- function(file_path,
                            demographic_col = "Black",
                            label_col = "Label",
                            sample_id_col = "sample_id_alias",
                            filter_label = "both",
                            demographic_group_in_bias_detection = c(1, -1),
                            demographic_mapping = c("-1" = "White", "1" = "Black", "0" = "White"),
                            label_mapping = c("-1" = "noncancer", "1" = "cancer")) {

  annot <- fread(file_path, header = TRUE, data.table = FALSE)
  colnames(annot)[str_detect(colnames(annot),"Label")] <- "Label"
 
  if (filter_label == "both") {
    annot <- annot[annot[[demographic_col]] %in% demographic_group_in_bias_detection & annot[[label_col]] == -1, ] ## allow demographic == 0
  } else if (filter_label == "demographic"){
    annot <- annot[annot[[demographic_col]] %in% demographic_group_in_bias_detection, ]
  } else if (filter_label == "none"){
    annot <- annot
  }
  
  annot[[demographic_col]] <- as.character(demographic_mapping[as.character(annot[[demographic_col]])])
  annot[[label_col]] <- as.character(label_mapping[as.character(annot[[label_col]])])
  
  annot$Group <- paste(annot[[demographic_col]], annot[[label_col]], sep = "-")
  annot
}

load_data <- function(data_file) {
  dt <- fread(data_file, header = TRUE, data.table = FALSE)
  colnames(dt)[1] <- "key"
  rownames(dt) <- dt$key
  dt
}

load_markers <- function(marker_file, dt) {
  markers <- fread(marker_file, header = FALSE, data.table = FALSE)$V1
  markers <- markers[markers %in% dt$key]
  return(markers)
}

