# ConsensusTranscriptomicSubtype

Classify patient samples into **Consensus Transcriptomic Subtypes (CTS)** of **sepsis** using packaged reference data, user-supplied healthy controls, **ComBat** normalization, and **random forest** classification.

The current implementation supports four CTS groups:

- **CTS1**
- **CTS2**
- **CTS3**
- **CTS4** — treated here as a **health/null subtype**

User-provided healthy controls are incorporated during normalization to preserve healthy-versus-sepsis biology while reducing dataset-specific technical effects.

---

## Highlights

- Packaged reference expression matrix and subtype labels
- Joint normalization of:
  - packaged reference samples
  - user healthy controls
  - user case samples
- CTS prediction with random forest classification
- Optional heatmap generation
- Optional silhouette calculation
- Optional PDF export of heatmap and silhouette plots

---

## Installation

### From source archive

```r
install.packages("ConsensusTranscriptomicSubtype_0.1.4.tar.gz", repos = NULL, type = "source")
````

If the archive is in another folder:

```r
install.packages(
  "C:/path/to/ConsensusTranscriptomicSubtype_0.1.4.tar.gz",
  repos = NULL,
  type = "source"
)
```

### From GitHub

Use this only if the repository contains the full package source files.

```r
# install.packages("remotes")
remotes::install_github("bpscicluna/ConsensusTranscriptomicSubtype")
```

---

## Load the package

```r
library(ConsensusTranscriptomicSubtype)
```

---

## Packaged reference data

The package includes:

* `exp_core_g` — reference gene expression matrix
* `core_samples` — reference sample annotations containing CTS labels

These datasets are loaded automatically by the classifier.

---

## Input requirements

The classifier expects **two user-provided datasets**:

* `new_case_data` — gene expression data for the samples to classify
* `new_healthy_data` — gene expression data for healthy control samples

### Accepted input types

Input can be provided as:

* `.csv`
* `.txt`
* `.tsv`
* `data.frame`
* `matrix`

### Required structure

For both `new_case_data` and `new_healthy_data`:

* genes in **rows**
* samples in **columns**
* expression values must be numeric
* gene identifiers supplied either:

  * in the **first column**, or
  * as **row names**

**Preferred identifier format:** Ensembl gene IDs

### Example input table

| gene       | sample_1 | sample_2 | sample_3 |
| ---------- | -------- | -------- | -------- |
| ENSG000001 | 5.3      | 4.8      | 6.1      |
| ENSG000002 | 2.1      | 2.5      | 1.9      |
| ENSG000003 | 8.0      | 7.6      | 7.8      |

---

## Read input files

### CSV

```r
new_case_data <- read.csv("your_case_samples.csv", check.names = FALSE)
new_healthy_data <- read.csv("your_healthy_controls.csv", check.names = FALSE)
```

### Tab-delimited

```r
new_case_data <- read.delim("your_case_samples.tsv", check.names = FALSE)
new_healthy_data <- read.delim("your_healthy_controls.tsv", check.names = FALSE)
```

---

## Main function

```r
run_subtype_classifier(
  new_case_data,
  new_healthy_data,
  gene_list = NULL,
  make_heatmap = TRUE,
  ntrees = 500,
  save_plots = FALSE,
  output_dir = ".",
  heatmap_file = "CTS_heatmap.pdf",
  silhouette_file = "CTS_silhouette.pdf"
)
```

### Arguments

| Argument           | Description                                                      |
| ------------------ | ---------------------------------------------------------------- |
| `new_case_data`    | Matrix or data frame containing the new case samples to classify |
| `new_healthy_data` | Matrix or data frame containing new healthy control samples      |
| `gene_list`        | Optional vector of genes to restrict the classifier              |
| `make_heatmap`     | Logical indicating whether to generate a heatmap                 |
| `ntrees`           | Number of trees for the random forest model                      |
| `save_plots`       | Logical indicating whether to save plots as PDF files            |
| `output_dir`       | Directory where plot files will be saved                         |
| `heatmap_file`     | File name for the heatmap PDF                                    |
| `silhouette_file`  | File name for the silhouette PDF                                 |

---

## Quick start

```r
library(ConsensusTranscriptomicSubtype)

new_case_data <- read.csv("your_case_samples.csv", check.names = FALSE)
new_healthy_data <- read.csv("your_healthy_controls.csv", check.names = FALSE)

result <- run_subtype_classifier(
  new_case_data = new_case_data,
  new_healthy_data = new_healthy_data
)

head(result$predictions)
```

---

## Save heatmap and silhouette plot

```r
result <- run_subtype_classifier(
  new_case_data = new_case_data,
  new_healthy_data = new_healthy_data,
  make_heatmap = TRUE,
  save_plots = TRUE,
  output_dir = "CTS_results",
  heatmap_file = "cts_heatmap.pdf",
  silhouette_file = "cts_silhouette.pdf"
)

result$plot_files
```

---

## Additional examples

### Disable heatmap generation

```r
result <- run_subtype_classifier(
  new_case_data = new_case_data,
  new_healthy_data = new_healthy_data,
  make_heatmap = FALSE
)
```

### Use a custom gene set

```r
selected_genes <- c("ENSG000001", "ENSG000002", "ENSG000003")

result <- run_subtype_classifier(
  new_case_data = new_case_data,
  new_healthy_data = new_healthy_data,
  gene_list = selected_genes
)
```

---

## Output

The function returns a list containing:

| Element                | Description                                                              |
| ---------------------- | ------------------------------------------------------------------------ |
| `predictions`          | Predicted CTS labels for the new case samples                            |
| `expression_corrected` | Corrected `reference`, `new_healthy`, and `new_case` matrices            |
| `silhouette`           | Silhouette object for the predicted case samples, when applicable        |
| `rf_model`             | Trained random forest model                                              |
| `genes_used`           | Genes used in classification                                             |
| `gene_overlap`         | Number of overlapping genes across reference, healthy, and case datasets |
| `batch`                | Batch labels used during ComBat normalization                            |
| `bio_group`            | Biological group labels used during ComBat normalization                 |
| `plot_files`           | File paths of exported PDF plots, if saved                               |

### Example output access

```r
result$predictions
result$gene_overlap
result$plot_files
```

---

## Example workflow

```r
library(ConsensusTranscriptomicSubtype)

# Read user data
new_case_data <- read.csv("your_case_samples.csv", check.names = FALSE)
new_healthy_data <- read.csv("your_healthy_controls.csv", check.names = FALSE)

# Run classifier and save visual output
result <- run_subtype_classifier(
  new_case_data = new_case_data,
  new_healthy_data = new_healthy_data,
  make_heatmap = TRUE,
  save_plots = TRUE,
  output_dir = "CTS_results",
  heatmap_file = "cts_heatmap.pdf",
  silhouette_file = "cts_silhouette.pdf"
)

# View predicted classes
result$predictions

# View saved plot locations
result$plot_files
```

---

## Interpretation notes

* **CTS4** is treated here as a **health/null subtype**
* user healthy controls are incorporated during normalization
* subtype classification is performed on the **new case samples**
* healthy controls should ideally come from the same study, platform, or preprocessing workflow as the case samples

---

## Best practices

For best performance:

* use healthy controls from the same study or platform as the case samples
* ensure gene identifiers match those used in the packaged reference data
* use Ensembl gene IDs where possible
* ensure sufficient overlap across:

  * packaged reference data
  * user healthy controls
  * user case samples

If gene overlap is too low, the function will stop with an error.

---

## Troubleshooting

### Too few overlapping genes across reference, new cases, and new healthy controls

Check that:

* gene IDs were imported correctly
* the first column containing gene IDs was not lost
* Ensembl IDs are used consistently

### `new_case_data` or `new_healthy_data` must have row names with gene IDs

Make sure gene identifiers are supplied either:

* in the first column of the file, or
* as row names after import

### No silhouette plot was saved

A silhouette plot is only produced when:

* there are enough classified case samples
* more than one predicted subtype is present

### Installation from `.tar.gz` fails on Windows

You may need the appropriate R build tools installed.

---

## Contributing

Bug reports, suggestions, and pull requests are welcome.
Please open an issue or pull request on GitHub.

---

## License

This package is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

See the full license text in the [LICENSE](LICENSE) file.

```
