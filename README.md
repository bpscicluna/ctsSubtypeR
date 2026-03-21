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

## Installation

### Option 1: install from the packaged archive (`.tar.gz`)

If you downloaded or uploaded the package as a source archive, install it in R using:

```r
install.packages("ctsSubtypeR_0.0.0.9000.tar.gz", repos = NULL, type = "source")


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

