# CrossBatch-ML

**Batch-aware validation for molecular classification under cohort and platform shift.**

## 🧠 Motivation

Random train/test splits can mix samples from different batches, cohorts, or platforms across both sets and therefore produce overly optimistic performance estimates.

The broader SSV study found a consistent pattern across multiple datasets and machine-learning methods:

> **As batch separation between training and testing increases, classification error tends to increase.**

This repository demonstrates that idea with one cross-platform expression dataset example.

---

## Validation Framework

The analysis compares progressively stronger batch-aware splits:

| Split | Interpretation |
|---|---|
| `50` | strongly mixed / random-like |
| `70` | partially separated |
| `80` | more separated |
| `90` | near-external validation |
| `100` | fully separated external validation |

Preprocessing is fitted on the training data and then applied to the test data to avoid train-test leakage.

The workflow evaluates:

- KNN
- LASSO
- PAM
- SVM
- Random Forest
- XGBoost

---

## Case Study

This repository uses datasets from:

**Gao et al. (2019), _DeepCC: a novel deep learning-based framework for cancer molecular subtype classification_.**  
https://doi.org/10.1038/s41389-019-0157-8

The original four-class CMS task is reformulated as:

```text
CMS4     → positive class
CMS1/2/3 → negative class
```

The source cohort is **TCGA RNA-seq**, and three independent GEO microarray cohorts are evaluated separately:

- GSE13067
- GSE13294
- GSE14333

For each TCGA–GEO comparison:

```text
gene harmonization
        ↓
TCGA-only feature selection
        ↓
sample-wise rank-normal alignment
        ↓
split = 50 / 70 / 80 / 90 / 100
        ↓
normalization × classifier evaluation
        ↓
classification error
```

The strictest setting (`100`) corresponds to **TCGA-only training and GEO-only testing**.

---

## Main Finding

Across all three external cohorts, classification error increased as TCGA and GEO became more completely separated between training and testing.

The deterioration was:

- moderate for GSE13067
- more consistent for GSE13294
- strongest for GSE14333

The alignment diagnostics show that obvious cross-platform scale differences are substantially reduced before classification, but stricter external validation still leads to higher error.

---

## Repository Structure

```text
CrossBatch-ML/
├── code/      # R workflows and data notes
├── figure/    # diagnostic and error-curve figures
├── README.md
└── .gitignore
```

Large expression matrices are not included in the repository. Data sources, expected file names, and preprocessing notes are documented in the README files under `code/train_s1/`, `code/validation_s1/`, and `code/label/`.

---

## Scope

This repository **does not reproduce the original DeepCC neural-network model**.

Instead, it uses datasets from that study as a real cross-platform example for evaluating how machine-learning validation changes under progressively stronger batch, cohort, and platform separation.
