# Introduction 
Analysis code (in R) for analysis of bulk and single-cell RNAseq data from pediatric leukemia patient samples.

In this study, we generated bulk and single-cell RNAseq profiles for a total of 69 pediatric leukemia patients treated at Children's Mercy Kansas City. The code here details data processing and differential expression analysis, as well as generating the figures for the manuscript (link to be generated). 

Manuscript available on [bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.23.696262v1).

# Data availability
Raw data (both bulk and single-cell) are available through [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002529.v2.p1). Processed single-cell data are available through the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP3385/single-cell-profiling-of-pediatric-leukemia-patient-samples).

All manuscript analysis and figures can be generated from the Seurat object downloadable through the SCP and input data provided in this repository.

# Inputs
The following files are in the inputs directory: 

- sampleTable.csv: sample metadata, also Supplementary Table 1 of the manuscript
- bulkRNA_kallisto_rawCounts.csv: raw counts for all bulk RNAseq samples
- pseudobulk_sampleCounts.csv: gene expression counts from scRNAseq, aggregated at the patient level for leukemic blasts only; code for generating this file from the Seurat object is included, but we are providing the counts for convenient further analysis

The Seurat object downloadable from the SCP corresponds to the 'merge_all.rds' object saved at the end of full_SO_setup.R, with minor changes to metadata made to accomodate SCP metadata conventions.

# Analysis files

Analysis was performed in the following order:

1. bulkRNA_leukemias.R - batch correction, visualization, and differential expression analysis across all bulk RNAseq samples
2. 10x outputs were processed using cellranger for each single-cell batch, read in, merged across pools, and filtered based on demuxlet singlet calls (see manuscript Materials & Methods)
3. full_SO_setup.R - reading in single-cell batch "raw" Seurat objects, cell type assignment, merging, dimensionality reduction
4. cancer_SO_setup.R - filtering full Seurat object down to presumed leukemic blasts only for downstream analysis; generation of pseudobulk counts for differential expression
5. scRNA_leukemias.R - pseudobulk differential expression and visualization across leukemias (scRNAseq)
6. B-ALL_subtypes.R - filtering down to B-ALLs only, differential expression across subtypes in bulk and sc RNAseq data
7. sc_overlapDE_plots.R - plotting differentially expressed genes (bulk and sc) across leukemias and B-ALL subtypes in the single-cell dataset
