# ============================================
# Command-line wrapper for DeBias pipeline
# ============================================

options(warn = -1)
suppressPackageStartupMessages(library(argparse))

# Source the actual pipeline function
source("07_run_pipeline.R")

# Define argument parser
parser <- ArgumentParser(description = "Run DeBias bias-removal pipeline")

parser$add_argument("--marker",         required = TRUE,  help = "Marker file path")
parser$add_argument("--annot_file",     required = TRUE,  help = "Annotation file path")
parser$add_argument("--data_file",      required = TRUE,  help = "Input data file path")
parser$add_argument("--output_dir",     required = TRUE,  help = "Output path")
parser$add_argument("--select_pcs_criterion",type = "character",default = "es",help = "Criterion for selecting PCs (e.g. es)")
parser$add_argument("--padj_cutoff",    type = "double",   default = 0.05,   help = "Adjusted p-value cutoff")
parser$add_argument("--pblack_cutoff",  type = "double",   default = 0.05,   help = "p-value cutoff for Black group")
parser$add_argument("--lda_p_cutoff",   type = "double",   default = 0.05,   help = "LDA p-value cutoff")
parser$add_argument("--es_race_cutoff", type = "double",   default = 0.1,    help = "Effect size cutoff for race")
parser$add_argument("--es_black_cutoff",type = "double",   default = 0.5,    help = "Effect size cutoff for Black group")
parser$add_argument("--lda_es_cutoff",  type = "double",   default = 0.1,    help = "LDA effect size cutoff")
parser$add_argument("--signal_test_set",type = "character",default = "full", help = "Testing set identifier")

# Parse arguments
args <- parser$parse_args()

# Call the main function
run_pipeline(
  marker_file          = args$marker,
  annot_file           = args$annot_file,
  data_file            = args$data_file,
  output_dir           = args$output_dir,
  select_pcs_criterion = args$select_pcs_criterion,
  padj_cutoff          = args$padj_cutoff,
  pblack_cutoff        = args$pblack_cutoff,
  lda_p_cutoff         = args$lda_p_cutoff,
  es_race_cutoff       = args$es_race_cutoff,
  es_black_cutoff      = args$es_black_cutoff,
  lda_es_cutoff        = args$lda_es_cutoff,
  signal_test_set      = args$signal_test_set
)


