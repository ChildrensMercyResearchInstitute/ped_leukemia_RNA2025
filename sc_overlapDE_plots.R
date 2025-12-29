library(tidyverse)
library(pheatmap)
library(Seurat)
library(stringr)
library(Cairo)
library(SingleCellExperiment)
library(Matrix.utils)
library(scCustomize)

options(future.globals.maxSize = 5e+10)

cancer_so <- readRDS('SeuratObjects/cancer_SO.rds')
cancer_so
#An object of class Seurat 
#69230 features across 124757 samples within 2 assays 
#Active assay: SCT (32629 features, 4304 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

#read in genes DE across bulk and scRNAseq, overlapping w/PMTL
goi <- read.csv('bulk_sc_PMTL_overlapGenes.csv', row.names = 1)

DefaultAssay(cancer_so) <- 'RNA'
Idents(cancer_so) <- cancer_so$leukemia

pdf('cancerSO_images/DotPlot_bulk_sc_PMTL_leukemias.pdf', height = 5, width = 12)
Clustered_DotPlot(cancer_so, features = rownames(goi), flip = TRUE, 
                  colors_use_exp = c('lightgrey', 'blue')) & RotatedAxis() 
dev.off()


### filter down to B-ALLs only -----

Idents(cancer_so) <- cancer_so$leukemia

BALL_so <- subset(cancer_so, idents = 'B-ALL')
BALL_so
#An object of class Seurat 
#69230 features across 85863 samples within 2 assays 
#Active assay: SCT (32629 features, 4304 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

DefaultAssay(BALL_so) <- 'RNA'
Idents(BALL_so) <- BALL_so$subtype

#remove NOSs (since not a unified class)
BALL_known_so <- subset(BALL_so, idents = '', invert = TRUE)
BALL_known_so
#An object of class Seurat 
#69230 features across 83989 samples within 2 assays 
#Active assay: RNA (36601 features, 0 variable features)
#2 layers present: counts, data
#1 other assay present: SCT
#2 dimensional reductions calculated: pca, umap

subtype_goi <- read.csv('subtype_bulk_sc_PMTL_overlapGenes.csv', row.names = 1)

pdf('cancerSO_images/B-ALL_DotPlot_bulk_sc_PMTL_leukemias.pdf', height = 5, width = 10)
Clustered_DotPlot(BALL_known_so, features = rownames(subtype_goi), flip = TRUE, 
                  colors_use_exp = c('lightgrey', 'blue')) & RotatedAxis() 
dev.off()

