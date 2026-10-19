# Validation Data

> **Note:** The validation expression files are too large to be included in this repository.  
> The original study cited below provides access to the corresponding datasets and download resources.

## Files

- `gse13067_rma_expression.tsv`
- `gse13294_rma_expression.tsv`
- `gse14333_rma_expression.tsv`

## Source

These datasets correspond to three external cohorts used in:

**Gao et al. (2019), _DeepCC: a novel deep learning-based framework for cancer molecular subtype classification_**  
[Oncogenesis, 8:44](https://doi.org/10.1038/s41389-019-0157-8)

The original datasets are GEO studies:

- GSE13067
- GSE13294
- GSE14333

In the DeepCC study, these cohorts were part of the curated colorectal cancer datasets obtained from the **CRCSC Synapse repository**:

[CRCSC Synapse Repository](https://www.synapse.org/#!Synapse:syn2623706/wiki/)

## Data Type

- Platform type: **microarray gene-expression data**
- Repository files used here: **RMA-processed expression matrices**
- Role in this project: **independent target / validation cohorts**

Each cohort is analyzed separately against the TCGA RNA-seq training cohort.

## Preprocessing

The files in this repository are treated as already processed and log-scale expression matrices.

In the current R workflow:

1. No additional log transformation is applied to the RMA validation matrices.
2. Microarray probe IDs are mapped to gene symbols.
3. Multiple probes mapping to the same gene are collapsed by mean expression.
4. Only genes shared with the TCGA training matrix are retained.
5. Cross-platform alignment is subsequently performed before classification.

The three validation cohorts follow the same preprocessing and evaluation workflow.
