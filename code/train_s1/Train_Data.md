# Training Data

This folder contains the source training expression dataset used in the R workflow.

## File

- `TCGACRC.txt` — TCGA colorectal cancer gene-expression data

## Source

The dataset follows the colorectal cancer analysis described in:

**Gao et al. (2019), _DeepCC: a novel deep learning-based framework for cancer molecular subtype classification_**  
[Oncogenesis, 8:44](https://doi.org/10.1038/s41389-019-0157-8)

The original TCGA expression data were obtained from the **Broad Firehose GDAC portal**.

## Data Type

- Platform: **TCGA RNA-seq**
- Expression type: gene-level RNA-seq expression
- Role in this project: **source / training cohort**

## Preprocessing

According to the original study:

1. Level-3 RNA-seq expression data were used.
2. RSEM scaled estimates were converted to TPM.
3. Expression values were log2-transformed.

In the current R workflow, `TCGACRC.txt` is read as the training expression matrix. If the values do not appear to be on a log scale, `log2(x + 1)` is applied before downstream processing.

Gene identifiers are subsequently harmonized with the validation datasets before classification.
