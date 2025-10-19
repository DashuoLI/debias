

# DeBias Pipeline

This repository contains a modular R pipeline for **DeBias**.

We present *DeBias*, a computational framework for mitigating demographic biases in high-dimensional biomedical datasets. *DeBias* identifies and removes bias-associated subspaces from the feature space using control samples, enabling global correction of demographic distortions while preserving disease-specific signals. To evaluate its effectiveness, we applied *DeBias* to cell-free DNA methylation data for cancer detection. *DeBias* achieved a significant reduction in the number of features exhibiting demographic bias and outperformed existing methods in improving cancer detection performance for minority populations. Performance gains were validated in independent cohorts, highlighting the robustness of the approach.

---

## 📂 Project Structure




├── 01_load_data.R        # Functions for loading data and annotations \
├── 02_preprocess.R       # Functions for preprocessing and filtering \
├── 03_pca_analysis.R     # PCA computation and biased PC selection \
├── 04_bias_model.R       # Iterative LDA for bias removal \
├── 05_apply_model.R      # Adjust rotation and reconstruct the adjusted matrix \
├── 06_postprocess.R      # Postprocessing and evaluation functions \
├── 07_run_pipeline.R     # Main script to run the entire pipeline \
├── data/                 # Example input data files \
├── output/               # Pipeline outputs \
└── README.md



---

## ⚙️ Installation

1. Clone the repository:

```bash
git clone https://github.com/yourusername/debias.git
cd debias
````

2. Install required R packages:

```r
install.packages(c("data.table", "MASS", "infotheo", "stringr", "coin", "effectsize"))
```

---

## 📝 Usage

### Running the pipeline from command line

The main script `debias.R` orchestrates all steps.

```bash
Rscript debias.R \
  --marker <marker_file> \
  --marker_type <marker_type> \
  --annot <annotation_file> \
  --split <split_name> \
  --data <data_file> \
  --select_pcs_criterion <p|es> \
  --padj_cutoff <numeric> \
  --pblack_cutoff <numeric> \
  --lda_p_cutoff <numeric> \
  --es_race_cutoff <numeric> \
  --es_black_cutoff <numeric> \
  --lda_es_cutoff <numeric> \
  --signal_test_set <minor|full>
```

**Arguments:**

| Argument                 | Description                                          |
| ------------------------ | ---------------------------------------------------- |
| `--marker`               | Marker file path (CSV)                               |
| `--marker_type`          | Type of markers (e.g., gene, protein)                |
| `--annot`                | Sample annotation file (CSV)                         |
| `--split`                | Split identifier (for multiple runs)                 |
| `--data`                 | Feature matrix file (CSV/RDS)                        |
| `--select_pcs_criterion` | Criterion for selecting biased PCs (`p` or `es`)     |
| `--padj_cutoff`          | Adjusted p-value cutoff for PC selection             |
| `--pblack_cutoff`        | p-value threshold for cancer signal testing          |
| `--lda_p_cutoff`         | LDA significance cutoff for bias removal             |
| `--es_race_cutoff`       | Effect size cutoff for bias detection                |
| `--es_black_cutoff`      | Effect size cutoff for cancer signal preservation    |
| `--lda_es_cutoff`        | LDA effect size cutoff for bias removal              |
| `--signal_test_set`      | Subset to use for signal testing (`minor` or `full`) |

---

## 🔹 Workflow Overview

1. **Load data and annotations** (`01_load_data.R`)
2. **Filter and preprocess** features (`02_preprocess.R`)
3. **Perform PCA**, identify biased PCs (`03_pca_analysis.R`)
4. **Iterative LDA** to remove demographic bias (`04_bias_model.R`)
5. **Adjust rotation and reconstruct** the corrected matrix (`05_apply_model.R`)
6. **Evaluate bias removal and signal preservation** (`06_postprocess.R`)
7. **Main execution** (`07_run_pipeline.R`)

---

## 📊 Output

The pipeline generates:

* **Adjusted feature matrix** (`output/adjusted_<marker_type>_<split>.csv.gz`)
* **Evaluation statistics** (`output/adjusted_<marker_type>_<split>.stats.txt`)
* Optional LDA/PC stats for downstream analysis


---

## 📜 Citation

If you use this pipeline in your research, please cite:

Shuo Li, Weihua Zeng, Wenyuan Li, Chun-Chi Liu, Yonggang Zhou, Xiaohui Ni, Mary L. Stackpole, Angela H. Yeh, Andrew Melehy, David S Lu, Steven S Raman, William Hsu, Lopa Mishra, Kirti Shetty, Benjamin Tran, Megumi Yokomizo, Preeti Ahuja, Yazhen Zhu, Hsian-Rong Tseng, Denise R. Aberle, Vatche G. Agopian, Steven-Huy B. Han, Samuel W. French, Steven M. Dubinett, Xianghong Jasmine Zhou, Wing Hung Wong. "Reducing demographic biases in biomedical machine learning: application to cancer detection using cfDNA methylation" Under Revision. 2025.



