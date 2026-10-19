# Code

This folder contains the R workflow used to evaluate classifier robustness under progressively stronger cohort and platform separation.

The analysis is based on the colorectal cancer case study from:

**Gao et al. (2019), _DeepCC: a novel deep learning-based framework for cancer molecular subtype classification_**  
[Oncogenesis, 8:44](https://doi.org/10.1038/s41389-019-0157-8)

The original four-class subtype problem was reformulated here as a binary classification task:

- `1` — CMS4
- `0` — Others (`CMS1 + CMS2 + CMS3`)

The main objective is to test how classification error changes when validation moves from partially mixed TCGA/GEO splits to strict external validation.

---

## Research Design

TCGA RNA-seq samples are treated as the source batch, while each GEO microarray cohort is treated as an independent target batch.

Three external validation analyses are performed separately:

- TCGA → GSE13067
- TCGA → GSE13294
- TCGA → GSE14333

These are parallel analyses rather than sequential validation steps.

The batch variable therefore represents both cohort origin and measurement platform.

---

## Workflow

For each TCGA–GEO comparison, the R workflow follows the same sequence:

1. Read TCGA RNA-seq, GEO microarray, and CMS label files.
2. Match sample identifiers with CMS labels.
3. Retain valid CMS1–CMS4 samples.
4. Convert the multiclass CMS labels to CMS4 vs Others.
5. Harmonize gene identifiers across TCGA and GEO.
6. Restrict both datasets to their common genes.
7. Select the top 1,500 genes by mean expression using TCGA only.
8. Perform sample-wise rank-normal alignment independently within each sample.
9. Remove invalid or zero-variance features.
10. Construct progressively stronger TCGA–GEO train/test separation.
11. Apply multiple downstream normalization methods.
12. Train and evaluate multiple machine-learning classifiers.
13. Summarize classification error across repeated splits.

The final aligned expression matrix is passed directly into `ml_sum()`.

---

## Train-Test Separation

The main analyses use:

```text
split = 50, 70, 80, 90, 100
