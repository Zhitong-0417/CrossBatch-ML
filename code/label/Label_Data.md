# CMS Labels

This folder contains the sample-level class labels used in the R workflow.

## File

- `cms_labels_public_all.txt`

## Source

The CMS classification labels follow the colorectal cancer datasets used in:

**Gao et al. (2019), _DeepCC: a novel deep learning-based framework for cancer molecular subtype classification_**  
[Oncogenesis, 8:44](https://doi.org/10.1038/s41389-019-0157-8)

The sample-level CMS labels were obtained from the **CRCSC Synapse repository**:

[CRCSC Synapse Repository](https://www.synapse.org/#!Synapse:syn2623706/wiki/)

## Data Type

The file contains sample identifiers, dataset identifiers, and consensus molecular subtype (CMS) assignments.

The original labels include:

- CMS1
- CMS2
- CMS3
- CMS4

## Use in This Project

For the current binary classification analysis, the original CMS labels are converted to:

- `1` — CMS4
- `0` — Others (`CMS1 + CMS2 + CMS3`)

Only samples with valid CMS1–CMS4 assignments are retained.

Sample identifiers are matched separately to the TCGA training cohort and each GEO validation cohort before model evaluation.
