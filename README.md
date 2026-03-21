# ctsSubtypeR

`ctsSubtypeR` is an R package for classifying new patient gene expression profiles into consensus transcriptomic subtypes (CTS) using packaged reference data, batch correction, and random forest classification.

The package currently supports four CTS groups:

- **CTS1**
- **CTS2**
- **CTS3**
- **CTS4**, defined here as a **health/null subtype**

Including CTS4 broadens the application of the framework beyond clearly diseased cohorts and may support future use in settings such as emergency medicine, screening, or primary care.

---

## What the package does

`ctsSubtypeR`:

- uses packaged reference expression data and subtype labels
- accepts new expression data as a matrix or data frame
- identifies overlapping genes between reference and new samples
- performs joint batch correction using `ComBat`
- trains a random forest classifier on the packaged reference cohort
- predicts CTS labels for new samples
- optionally generates a heatmap of classified samples
- returns corrected matrices and classification outputs for downstream analysis

---

## Package contents

The package includes:

- **`run_subtype_classifier()`**  
  Main function for classifying new samples

- **`exp_core_g`**  
  Packaged reference gene expression matrix

- **`core_samples`**  
  Packaged reference sample annotation data frame containing the `CTS` labels

These packaged datasets are used automatically by the classifier and do not need to be provided by the user.

---

How to use the package
Basic usage
library(ctsSubtypeR)

new_expr <- read.csv("your_new_data.csv", check.names = FALSE)

result <- run_subtype_classifier(
  new_expr_data = new_expr
)

head(result$predictions)
Disable the heatmap
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  make_heatmap = FALSE
)
Use a custom number of trees
result <- run_subtype_classifier(
  new_expr_data = new_expr,
  ntrees = 1000
)
Restrict the classifier to a selected gene set
selected_genes <- c("ENSG000001", "ENSG000002", "ENSG000003")

result <- run_subtype_classifier(
  new_expr_data = new_expr,
  gene_list = selected_genes
)
Output

The function returns a list with the following components.

predictions

A data frame containing the predicted subtype for each sample.

result$predictions

Typical structure:

sample_id	CTS
sample_1	1
sample_2	4
sample_3	2
expression_corrected

A list containing the batch-corrected expression matrices:

core = corrected packaged reference matrix
new_data = corrected new-sample matrix
result$expression_corrected$core
result$expression_corrected$new_data
silhouette

A silhouette object describing how consistently the predicted new samples group according to their assigned subtype.

result$silhouette
rf_model

The fitted random forest model trained on the packaged reference data.

result$rf_model
genes_used

A character vector containing the genes used in the classification after intersecting reference and new data.

result$genes_used
gene_overlap

The number of overlapping genes between the packaged reference dataset and the new input data.

result$gene_overlap
Example workflow
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
Interpretation of CTS labels

The package currently supports four CTS labels:

CTS1
CTS2
CTS3
CTS4

Within the current framework, CTS4 is treated as a health/null subtype. This does not necessarily imply absence of all biology, but rather a reference-like or non-disease-associated class within the current model.

Important input considerations

For best performance:

gene identifiers in the new data should match those used in the packaged reference data
gene identifiers should ideally be Ensembl IDs
row names should contain gene IDs whenever possible
sample IDs should be unique
gene IDs should be unique
the new dataset must contain sufficient overlap with the reference gene set

If too few overlapping genes are found, the function will stop and return an error.

Troubleshooting
“Too few overlapping genes between reference and new data”

Your input gene identifiers likely do not match the reference gene IDs used in the package.

Check:

whether you are using Ensembl IDs or gene symbols
whether row names were imported correctly
whether the first column containing genes was properly interpreted
“new_expr_data must have row names containing gene IDs”

Make sure gene identifiers are stored as row names, or provide them in the first column of the data frame.

Heatmap is cluttered

For larger datasets, try:

result <- run_subtype_classifier(
  new_expr_data = new_expr,
  make_heatmap = FALSE
)
Installation from .tar.gz fails on Windows

You may need the appropriate R build tools installed if installing from source.

Author

Brendon Scicluna

License

MIT


## Required input file type

The package expects a **text-based gene expression table** that can be read into R as either:

- a **CSV file** (`.csv`)
- a **tab-delimited text file** (`.txt` or `.tsv`)
- an **R matrix**
- an **R data frame**

For most users, the recommended input format is:

- **`.csv` file**
- genes in **rows**
- samples in **columns**
- first column containing **gene identifiers**

### Recommended file format

| gene | sample_1 | sample_2 | sample_3 |
|------|----------|----------|----------|
| ENSG000001 | 5.3 | 4.8 | 6.1 |
| ENSG000002 | 2.1 | 2.5 | 1.9 |
| ENSG000003 | 8.0 | 7.6 | 7.8 |

### Accepted gene identifier placement

Gene IDs can be supplied in either of these ways:

1. **Preferred:** as the first column of the file  
2. As **row names** after import into R

### Recommended gene identifier type

Gene identifiers should match the reference data used in the package.  
The safest option is to use:

- **Ensembl gene IDs**

### Reading input files into R

#### CSV file

```r
new_expr <- read.csv("your_new_data.csv", check.names = FALSE)

