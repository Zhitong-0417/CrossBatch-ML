# Code

This folder contains the R implementation of the cross-batch validation analysis used in this project.

## Research Objective

The goal of this analysis is to evaluate how classifier performance changes as the training and test sets become progressively more separated by cohort and measurement platform.

The original four-class subtype problem is reformulated as a binary classification task:

- **Positive class:** CMS4
- **Negative class:** CMS1 + CMS2 + CMS3

The source cohort is TCGA RNA-seq, and three independent GEO microarray cohorts are used as external validation targets:

- GSE13067
- GSE13294
- GSE14333

Each GEO cohort is analyzed independently against TCGA, so the three analyses are parallel cross-cohort comparisons rather than a sequential validation procedure.

The central question is whether validation performance becomes more pessimistic as the train/test split moves from batch-mixed validation toward fully separated external validation.

---

## Analysis Framework

For each TCGA–GEO comparison, the workflow follows the same structure:

1. **Load expression and CMS label data**
2. **Define the binary outcome**
   - CMS4 = 1
   - CMS1/CMS2/CMS3 = 0
3. **Define the batch/domain variable**
   - Batch 1 = TCGA RNA-seq
   - Batch 2 = one GEO microarray cohort
4. **Harmonize gene identifiers**
   - TCGA and GEO data are converted to gene-level matrices
   - Only genes shared by the two cohorts are retained
5. **Training-only feature selection**
   - Features are ranked using TCGA mean expression only
   - The top 1,500 genes are retained
6. **Cross-platform alignment**
   - TCGA and GEO samples are transformed independently using sample-wise rank-normal alignment
   - Within each sample, gene ranks are mapped to approximately standard-normal quantiles
7. **Remove invalid or zero-variance features**
8. **Construct progressively stronger batch-separated train/test splits**
9. **Evaluate multiple normalization methods and classifiers**
10. **Summarize classification error across repeated runs**

---

## Train/Test Split Design

The main SSV analysis uses:

```r
split_combination <- c(50, 70, 80, 90, 100)
```

With TCGA treated as the source batch and GEO as the target batch:

| Split | Training composition | Test composition | Interpretation |
|---|---|---|---|
| 50 | 50% of TCGA + 50% of GEO | remaining 50% of TCGA + 50% of GEO | strongly mixed |
| 70 | 70% of TCGA + 30% of GEO | 30% of TCGA + 70% of GEO | partially separated |
| 80 | 80% of TCGA + 20% of GEO | 20% of TCGA + 80% of GEO | more separated |
| 90 | 90% of TCGA + 10% of GEO | 10% of TCGA + 90% of GEO | near-external validation |
| 100 | TCGA only | GEO only | strict external validation |

Because the TCGA and GEO cohorts can have different sample sizes, these percentages refer to the fraction sampled **within each batch**, not necessarily the final proportion of the combined training set.

The strictest setting (`split = 100`) corresponds to:

```text
Train = TCGA only
Test  = GEO only
```

---

## Cross-Platform Alignment

Before classification, TCGA and GEO samples are aligned independently using a sample-wise rank-normal transformation.

For each sample:

\[
r_j = \operatorname{rank}(x_j)
\]

\[
p_j = \frac{r_j - 0.5}{m}
\]

\[
z_j = \Phi^{-1}(p_j)
\]

where \(m\) is the number of valid features in the sample.

This transformation places samples from different platforms on a common rank-based scale without using the test-set distribution to normalize the training data.

The `figure/` directory contains diagnostic plots before and after alignment.

---

## Downstream Normalization

After cross-platform alignment, the resulting continuous expression matrix is passed to the non-RNA-seq branch of `ml_sum()`.

The following downstream normalization settings are evaluated:

- `non` — no additional normalization
- `qn` — quantile normalization
- `mn` — median normalization
- `vsn` — variance stabilizing normalization

These transformations are fitted using the training data and then applied to the test data within the workflow.

---

## Classifiers

The following classifiers are evaluated:

- K-nearest neighbors (`knn`)
- Lasso logistic regression (`lasso`)
- Support vector machine (`svm`)
- Prediction Analysis for Microarrays (`pam`)
- Random forest (`rf`)
- XGBoost (`xgb`)

Each model is repeatedly evaluated across the predefined split settings.

---

## Main Scripts

### `SSV summary.R`

Core analysis utilities shared by all dataset-specific workflows.

Main functions include:

- normalization routines
- classifier training and prediction
- `ml_sum()` for repeated split-based evaluation
- `result_plot2()` for summarizing and plotting error curves

---

### `DeepCC_GSE13067_CMS4_ml_sum_strict.Rmd`

Main TCGA-to-GSE13067 analysis.

Implements:

```text
TCGA RNA-seq
    ↓
gene harmonization
    ↓
TCGA-only feature selection
    ↓
sample-wise rank-normal alignment
    ↓
split 50 / 70 / 80 / 90 / 100
    ↓
normalization × classifier evaluation
    ↓
classification error
```

---

### `DeepCC_GSE13294_CMS4_ml_sum_strict.Rmd`

Main TCGA-to-GSE13294 analysis using the same workflow.

The analysis is run independently from the GSE13067 comparison.

---

### `DeepCC_GSE14333_CMS4_ml_sum_strict.Rmd`

Main TCGA-to-GSE14333 analysis using the same workflow.

This cohort shows the largest deterioration in classification performance under strict TCGA-to-GEO separation among the three external cohorts.

---

### `DeepCC_GSE13067_GSE13294_CMS4_internal_platform_*.Rmd`

Supplementary analysis comparing internal platform/cohort structure.

This script is separate from the three main TCGA-to-GEO external validation analyses and is not required to reproduce the primary results.

---

## Diagnostic Plots

Before classification, sample-level expression diagnostics are generated to visualize platform-associated distribution differences.

### Before alignment

Low-to-middle-expression genes are used as poor-marker diagnostic features to visualize differences in scale, spread, and dynamic range between TCGA and GEO samples.

### After alignment

A separate diagnostic subset is selected from the aligned model-input features.

These plots are used to confirm that the sample-wise rank-normal transformation substantially reduces the most obvious cross-platform scale and spread differences.

The pre- and post-alignment plots do **not** use exactly the same genes and therefore should be interpreted as distribution diagnostics rather than a matched-feature before/after experiment.

---

## Main Result

Across GSE13067, GSE13294, and GSE14333, classification error generally increases as the train/test split becomes more strongly separated by cohort.

The main pattern is:

```text
batch-mixed validation
        ↓
increasing cohort/platform separation
        ↓
strict TCGA-only → GEO-only validation
        ↓
higher classification error
```

The deterioration is:

- moderate for GSE13067
- more consistent for GSE13294
- strongest for GSE14333

These results show that mixed-cohort splits can produce more optimistic error estimates than fully separated external validation.

---

## Data

The expression datasets are not stored in this folder.

See the README files under:

```text
train_s1/
validation_s1/
label/
```

for data sources, file names, preprocessing notes, and expected local file structure.

---

## Reference

Gao F, Wang W, Tan M, et al.  
**DeepCC: a novel deep learning-based framework for cancer molecular subtype classification.**  
*Oncogenesis*. 2019;8:44.  
https://doi.org/10.1038/s41389-019-0157-8
