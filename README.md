# ctsSubtypeR

`ctsSubtypeR` is an R package for classifying new patient gene expression profiles into consensus transcriptomic subtypes (CTS) using packaged reference data, batch correction, and random forest modelling.

The package currently supports four CTS groups:

- **CTS1**
- **CTS2**
- **CTS3**
- **CTS4**, defined here as a **health/null subtype**

The inclusion of CTS4 broadens the potential use of the framework beyond clearly diseased cohorts and may support future application in emergency medicine, screening, or primary care settings.

---

## Overview

`ctsSubtypeR` is designed to classify new transcriptomic profiles against an internal reference dataset packaged with the library. The workflow:

1. loads the packaged reference expression matrix and subtype labels
2. accepts user-supplied new expression data
3. matches overlapping genes between reference and new data
4. performs joint batch correction using `ComBat`
5. trains a random forest classifier on the reference cohort
6. predicts CTS labels for the new samples
7. optionally generates a heatmap of classified samples
8. returns subtype predictions, corrected matrices, and classification metrics

---

## Package contents

The package includes:

- **`run_subtype_classifier()`**  
  Main function for classifying new samples

- **`exp_core_g`**  
  Packaged reference gene expression matrix

- **`core_samples`**  
  Packaged reference sample annotation data frame containing the `CTS` labels

These packaged datasets are used automatically by the classifier and do not need to be supplied by the user.

---

## Installation

### Option 1: install from the packaged archive (`.tar.gz`)

If you downloaded or uploaded the package as a source archive, install it in R with:

```r
install.packages("ctsSubtypeR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
````

If the file is stored in another folder, provide the full path:

```r
install.packages("C:/path/to/ctsSubtypeR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

### Option 2: install from GitHub source

Use this only if the GitHub repository contains the full package source files such as `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, and `data/`.

```r
# install.packages("remotes")
remotes::install_github("bpscicluna/ctsSubtypeR")
```

---

## Loading the package

After installation, load the package in R:

```r
library(ctsSubtypeR)
```

---

## Required input file type

The package expects a gene expression table that can be read into R as either:

* a **CSV file** (`.csv`)
* a **tab-delimited text file** (`.txt` or `.tsv`)
* an **R matrix**
* an **R data frame**

For most users, the recommended format is a **CSV file** with:

* genes in **rows**
* samples in **columns**
* the first column containing **gene identifiers**

---

## Input data format

The function expects **gene expression data with genes in rows and samples in columns**.

### Required structure

* genes in **rows**
* samples in **columns**
* expression values must be numeric
* gene identifiers should be supplied either:

  * in the **first column** of the file, or
  * as **row names** after import into R

### Recommended gene identifier type

Gene identifiers should match those used in the packaged reference data.
The safest option is to use:

* **Ensembl gene IDs**

### Example input table

| gene       | sample_1 | sample_2 | sample_3 |
| ---------- | -------- | -------- | -------- |
| ENSG000001 | 5.3      | 4.8      | 6.1      |
| ENSG000002 | 2.1      | 2.5      | 1.9      |
| ENSG000003 | 8.0      | 7.6      | 7.8      |

---

## Reading the input file

### CSV file

```r
new_expr <- read.csv("your_new_data.csv", check.names = FALSE)
```

### Tab-delimited file

```r
new_expr <- read.delim("your_new_data.tsv", check.names = FALSE)
```

or

```r
new_expr <- read.delim("your_new_data.txt", check.names = FALSE)
```

---

## Main function

```r
run_subtype_classifier(
  new_expr_data,
  gene_list = NULL,
  make_heatmap = TRUE,
  ntrees = 500
)
```

### Arguments

* **`new_expr_data`**
  A matrix or data frame containing the new samples to classify

* **`gene_list`**
  Optional vector of gene identifiers used to restrict the classifier to a predefined gene set. If `NULL`, all genes in the packaged reference dataset are considered.

* **`make_heatmap`**
  Logical value indicating whether a heatmap of classified samples should be generated

* **`ntrees`**
  Number of trees used in the random forest classifier

---

## How to use the package

### Basic usage

```r
library(ctsSubtypeR)

new_expr <- read.csv("your_new_data.csv", check.names = FALSE)

result <- run_subtype_classifier(
  new_expr_data = new_expr
)

head(result$predictions)
```

### Disable the heatmap

```r
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  make_heatmap = FALSE
)
```

### Use a custom number of trees

```r
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  ntrees = 1000
)
```

### Restrict the classifier to a selected gene set

```r
selected_genes <- c("ENSG000001", "ENSG000002", "ENSG000003")

result <- run_subtype_classifier(
  new_expr_data = new_expr,
  gene_list = selected_genes
)
```

---

## Output

The function returns a list with the following components.

### `predictions`

A data frame containing the predicted subtype for each sample.

```r
result$predictions
```

Typical structure:

| sample_id | CTS |
| --------- | --- |
| sample_1  | 1   |
| sample_2  | 4   |
| sample_3  | 2   |

### `expression_corrected`

A list containing the batch-corrected expression matrices:

* `core` = corrected packaged reference matrix
* `new_data` = corrected new-sample matrix

```r
result$expression_corrected$core
result$expression_corrected$new_data
```

### `silhouette`

A silhouette object describing how consistently the predicted new samples group according to their assigned subtype.

```r
result$silhouette
```

### `rf_model`

The fitted random forest model trained on the packaged reference data.

```r
result$rf_model
```

### `genes_used`

A character vector containing the genes used in the classification after intersecting reference and new data.

```r
result$genes_used
```

### `gene_overlap`

The number of overlapping genes between the packaged reference dataset and the new input data.

```r
result$gene_overlap
```

---

## Example workflow

```r
library(ctsSubtypeR)

# Read new expression data
new_expr <- read.csv("your_new_data.csv", check.names = FALSE)

# Run classifier
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  make_heatmap = TRUE,
  ntrees = 500
)

# View subtype predictions
print(result$predictions)

# Check overlap with reference genes
print(result$gene_overlap)

# Access corrected expression matrix for new samples
head(result$expression_corrected$new_data)
```

---

## Interpretation of CTS labels

The package currently supports four CTS labels:

* **CTS1**
* **CTS2**
* **CTS3**
* **CTS4**

Within the current framework, **CTS4** is treated as a **health/null subtype**. This does not necessarily imply absence of all biology, but rather a reference-like or non-disease-associated class within the current model.

---

## Important input considerations

For best performance:

* gene identifiers in the new data should match those used in the packaged reference data
* gene identifiers should ideally be Ensembl IDs
* row names should contain gene IDs whenever possible
* sample IDs should be unique
* gene IDs should be unique
* the new dataset must contain sufficient overlap with the reference gene set

If too few overlapping genes are found, the function will stop and return an error.

---

## Troubleshooting

### “Too few overlapping genes between reference and new data”

Your input gene identifiers likely do not match the reference gene IDs used in the package.

Check:

* whether you are using Ensembl IDs or gene symbols
* whether row names were imported correctly
* whether the first column containing genes was properly interpreted

### “`new_expr_data` must have row names containing gene IDs”

Make sure gene identifiers are stored as row names, or provide them in the first column of the data frame.

### Heatmap is cluttered

For larger datasets, try:

```r
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  make_heatmap = FALSE
)
```

### Installation from `.tar.gz` fails on Windows

You may need the appropriate R build tools installed if installing from source.

---

## Author

Brendon Scicluna

---

## License

MIT

```

Save that as `README.md` in the root of your GitHub repository.
```
