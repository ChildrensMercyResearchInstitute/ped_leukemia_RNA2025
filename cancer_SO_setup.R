library(tidyverse)
library(pheatmap)
library(Seurat)
library(stringr)
library(Cairo)
library(SingleCellExperiment)
library(Matrix.utils)
library(scales)

options(future.globals.maxSize = 5e+10)

so <- readRDS('SeuratObjects/merge_all.rds')
so
#An object of class Seurat 
#69230 features across 144747 samples within 2 assays 
#Active assay: SCT (32629 features, 4304 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

### filter down to leukemic blasts only -----

Idents(so) <- so$seurat_clusters
# B-cells: 12, 14, 24, 40
# T-cells/NK cells: 17, 36
# Monocytes: 10, 34
# Erythrocytes: 29
# Macrophages: 32

cancer_so <- subset(so, idents = c(10, 12, 14, 17, 24, 29, 32, 34, 36, 40), invert = TRUE)
cancer_so
#An object of class Seurat 
#69230 features across 124757 samples within 2 assays 
#Active assay: SCT (32629 features, 4304 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

## plot metadata
to_plot <- c('batch', 'sample', 'patient', 'tissue', 'timepoint', 'leukemia', 
             'subtype', 'relapse', 'Phase', 'SingleR_hpca', 'SingleR_encode')

for(val in to_plot){
  name_nolab <- paste0('cancerSO_images/DimPlot_', val, '_nolab.pdf')
  cairo_pdf(name_nolab, width = 6, height = 6)
  print(DimPlot(cancer_so, group.by = val) & NoLegend())
  dev.off()
  name_lab <- paste0('cancerSO_images/DimPlot_', val, '.pdf')
  cairo_pdf(name_lab, width = 8, height = 6)
  print(DimPlot(cancer_so, group.by = val))
  dev.off()
}

## save Seurat object
saveRDS(cancer_so, 'SeuratObjects/cancer_SO.rds')


### Extract pseudobulk counts on a per-sample basis -----

# note that this is purposely pulling only counts for leukemic blasts
# for downstream pseudobulk differential expression analysis

counts <- cancer_so[["RNA"]]$counts
metadata <- cancer_so@meta.data

#create sce object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
sce
#class: SingleCellExperiment 
#dim: 36601 124757 
#metadata(0):
#  assays(1): counts
#rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
#rowData names(0):
#  colnames(124757): AAACCTGGTAGAGTGC-1_1 AAACCTGGTGCTAGCC-1_1 ... TTTGTTGGTTAGCATG-1-20
#TTTGTTGGTTATCCTA-1-20
#colData names(29): orig.ident nCount_RNA ... seurat_clusters SCT_snn_res.0.8
#reducedDimNames(0):
#  mainExpName: NULL
#altExpNames(0):

#want to aggregate by sample
samples <- purrr::set_names(levels(as.factor(sce$sample)))
n_cells <- as.numeric(table(sce$sample))  
m <- match(samples, sce$sample)  
ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL)  
dim(ei)  
#70 30

write.csv(ei, 'pseudobulk_cancer/pseudobulk_sampleTable.csv', row.names = FALSE, quote = FALSE)

#aggregate counts
groups <- colData(sce)[ , 'sample']
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum")
dim(pb)  
#70 36601

sample_counts <- t(pb)

write.csv(sample_counts, 'pseudobulk_cancer/pseudobulk_sampleCounts.csv', row.names = TRUE, quote = FALSE)


## leukemia DimPlot redo to use leukemia color scheme

#read in sampleTable
sampleTable <- read.csv('inputs/sampleTable.csv')

leukemia_color_list <- hue_pal()(7)
names(leukemia_color_list) <- unique(sampleTable$leukemia)

cairo_pdf('cancerSO_images/DimPlot_leukemia_nolab.pdf', width = 6, height = 6)
print(DimPlot(cancer_so, group.by = 'leukemia', cols = leukemia_color_list) & NoLegend())
dev.off()
cairo_pdf('cancerSO_images/DimPlot_leukemia.pdf', width = 8, height = 6)
print(DimPlot(cancer_so, group.by = 'leukemia', cols = leukemia_color_list))
dev.off()

## SingleR_encode DimPlot redo to use full color scheme

singler_color_list <- hue_pal()(16)
names(singler_color_list) <- unique(so$SingleR_encode)[order(unique(so$SingleR_encode))]

cairo_pdf('cancerSO_images/DimPlot_SingleR_encode_nolab.pdf', width = 6, height = 6)
print(DimPlot(cancer_so, group.by = 'SingleR_encode', cols = singler_color_list) & NoLegend())
dev.off()
cairo_pdf('cancerSO_images/DimPlot_SingleR_encode.pdf', width = 8, height = 6)
print(DimPlot(cancer_so, group.by = 'SingleR_encode', cols = singler_color_list))
dev.off()


