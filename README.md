# skin-sex

Code accompanying the manuscript:

**Lineage-structured sexual dimorphism in healthy human truncal skin revealed by single-cell transcriptomics**

## Overview

This repository contains the code used for the analysis of single-cell RNA-seq data from healthy adult human truncal skin, with a focus on sex-associated differences in cellular composition and lineage-specific transcriptional programs.

The main analyses include:

* preprocessing and quality control
* batch correction and integration
* donor-level cell composition analysis
* pseudo-bulk differential expression analysis
* Gene Ontology enrichment analysis
* cross-cell-type sharing analysis

```
code/
├── 00_cellranger.sh
├── 01_Data_import_and_Processing.py
├── 02_decontx.R
├── 03_Basic_Process_for_the_remove_ambientRNA_count.py
├── 04_integration_scvi.py
├── 05_composition_analysis.R
├── 06_pseudobulk_DE.R
├── 07_GO.R
└── 08_cross_celltype_sharing.R
```


## Requirements

Analyses were performed using:

* Cell Ranger: v6.1.2
* Python: 3.9 (Scanpy v1.9.3; scvi-tools v0.19.0; Scrublet v0.2.3)
* R: 4.2.1 (DaMiRseq v2.8.0; limma v3.52.4;  clusterProfiler v4.4.4; ComplexHeatmap v2.15.4)

Main packages include Scanpy, scvi-tools, Scrublet, limma, DaMiRseq, clusterProfiler, ComplexHeatmap, and ggplot2.

## Notes

Some scripts may require local path modification before execution. Please update file paths according to your computing environment.

## Contact

For questions regarding code or analysis, please contact the corresponding authors listed in the manuscript.
