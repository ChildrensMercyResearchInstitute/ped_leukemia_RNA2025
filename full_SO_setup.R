library(tidyverse)
library(pheatmap)
library(Seurat)
library(stringr)
library(Cairo)
library(scales)

library(SingleCellExperiment)
library(SingleR)
library(celldex)

options(future.globals.maxSize = 5e+10)

#read in sampleTable
sampleTable <- read.csv('inputs/sampleTable.csv')


### Process SO for pilot batch -----

pilot_so <- readRDS('pilot_merged_SeuratObj.rds')
pilot_so
#An object of class Seurat 
#36601 features across 32076 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#2 layers present: counts, data

Idents(pilot_so) <- pilot_so$orig.ident

pdf('pilot_images/VlnQC_raw.pdf', height = 4, width = 6)
VlnPlot(pilot_so, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), 
        pt.size = 0) & NoLegend()
dev.off()

## filter

pilot_so <- subset(pilot_so,
                   subset = nCount_RNA > 850 & 
                     nCount_RNA < 25000 &
                     percent_mt < 20)

pilot_so
#An object of class Seurat 
#36601 features across 31554 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#2 layers present: counts, data

pdf('pilot_images/VlnQC_filter.pdf', height = 4, width = 6)
VlnPlot(pilot_so, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), 
        pt.size = 0) & NoLegend()
dev.off()


## normalize RNA assay (for downstream visualization)
pilot_so <- NormalizeData(pilot_so)

## SCTransform (for dimensionality reduction, clustering)
pilot_so <- SCTransform(pilot_so)

## Dimensionality reduction
pilot_so <- RunPCA(pilot_so)
pilot_so <- RunUMAP(pilot_so, dims = 1:30)

## Add sample metadata
pilot_so$batch <- 'pilot'
pilot_so$patient <- substr(pilot_so$sample, 1, 9)
pilot_so$tissue <- sampleTable$tissue[match(pilot_so$sample, sampleTable$ID)]
pilot_so$timepoint <- sampleTable$timepoint[match(pilot_so$sample, sampleTable$ID)]
pilot_so$leukemia <- sampleTable$leukemia[match(pilot_so$sample, sampleTable$ID)]
pilot_so$subtype <- sampleTable$subtype[match(pilot_so$sample, sampleTable$ID)]
pilot_so$relapse <- sampleTable$relapse[match(pilot_so$sample, sampleTable$ID)]


## Cell cycle & cell type
DefaultAssay(pilot_so) <- "RNA"

#cell cycle data
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#cell type data (celldex)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
encode.se <- celldex::BlueprintEncodeData()

## assign cell cycle scores
pilot_so <- CellCycleScoring(pilot_so, s.features = s.genes, 
                             g2m.features = g2m.genes, set.ident = TRUE)
pilot_so$CC.difference <- pilot_so$S.Score - pilot_so$G2M.Score

## assign cell types
single_cell_experiment_object <- as.SingleCellExperiment(pilot_so)
predicted_cell_types_HPCA <- SingleR(single_cell_experiment_object, 
                                     ref = hpca.se, labels = hpca.se$label.main)
predicted_cell_types_ENCODE <- SingleR(single_cell_experiment_object,
                                       ref = encode.se, labels = encode.se$label.main)
pilot_so$SingleR_hpca <- predicted_cell_types_HPCA$labels
pilot_so$SingleR_encode <- predicted_cell_types_ENCODE$labels


## plot metadata
to_plot <- c('sample', 'patient', 'tissue', 'timepoint', 'leukemia', 
             'subtype', 'relapse', 'Phase', 'SingleR_hpca', 'SingleR_encode')

for(val in to_plot){
  name_nolab <- paste0('pilot_images/DimPlot_', val, '_nolab.pdf')
  cairo_pdf(name_nolab, width = 6, height = 6)
  print(DimPlot(pilot_so, group.by = val) & NoLegend())
  dev.off()
  name_lab <- paste0('pilot_images/DimPlot_', val, '.pdf')
  cairo_pdf(name_lab, width = 8, height = 6)
  print(DimPlot(pilot_so, group.by = val))
  dev.off()
}

## save Seurat object
saveRDS(pilot_so, 'SeuratObjects/pilot_processed.rds')

## pull out sample-specific barcodes for subset-bam
# NOTE: actually don't need these for pilot - already in CCDI

pilot_metadata <- pilot_so@meta.data

for(sample in unique(pilot_so$sample)){
  list_name <- paste0('pilot_barcodes/', sample, '.tsv')
  print(list_name)
  subset <- pilot_metadata[pilot_metadata$sample == sample, ]
  list <- substr(rownames(subset), 1, 18)
  print(length(list))
  write.table(list, list_name, sep = '\t', 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


### Process SO for s2022 batch -----

s2022_so <- readRDS('s2022_merged_SeuratObj.rds')
s2022_so
#An object of class Seurat 
#36601 features across 126402 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#2 layers present: counts, data

Idents(s2022_so) <- s2022_so$orig.ident

pdf('s2022_images/VlnQC_raw.pdf', height = 4, width = 6)
VlnPlot(s2022_so, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), 
        pt.size = 0) & NoLegend()
dev.off()

## filter

s2022_so <- subset(s2022_so,
                   subset = nCount_RNA > 850 & 
                     nCount_RNA < 25000 &
                     percent_mt < 20)

s2022_so
#An object of class Seurat 
#36601 features across 113193 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#2 layers present: counts, data

pdf('s2022_images/VlnQC_filter.pdf', height = 4, width = 6)
VlnPlot(s2022_so, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), 
        pt.size = 0) & NoLegend()
dev.off()


## normalize RNA assay (for downstream visualization)
s2022_so <- NormalizeData(s2022_so)

## SCTransform (for dimensionality reduction, clustering)
s2022_so <- SCTransform(s2022_so)

## Dimensionality reduction
s2022_so <- RunPCA(s2022_so)
s2022_so <- RunUMAP(s2022_so, dims = 1:30)

## Add sample metadata
s2022_so$batch <- 's2022'
s2022_so$patient <- substr(s2022_so$sample, 1, 9)
s2022_so$tissue <- sampleTable$tissue[match(s2022_so$sample, sampleTable$ID)]
s2022_so$timepoint <- sampleTable$timepoint[match(s2022_so$sample, sampleTable$ID)]
s2022_so$leukemia <- sampleTable$leukemia[match(s2022_so$sample, sampleTable$ID)]
s2022_so$subtype <- sampleTable$subtype[match(s2022_so$sample, sampleTable$ID)]
s2022_so$relapse <- sampleTable$relapse[match(s2022_so$sample, sampleTable$ID)]


## Cell cycle & cell type
DefaultAssay(s2022_so) <- "RNA"

## assign cell cycle scores
s2022_so <- CellCycleScoring(s2022_so, s.features = s.genes, 
                             g2m.features = g2m.genes, set.ident = TRUE)
s2022_so$CC.difference <- s2022_so$S.Score - s2022_so$G2M.Score

## assign cell types
single_cell_experiment_object <- as.SingleCellExperiment(s2022_so)
predicted_cell_types_HPCA <- SingleR(single_cell_experiment_object, 
                                     ref = hpca.se, labels = hpca.se$label.main)
predicted_cell_types_ENCODE <- SingleR(single_cell_experiment_object,
                                       ref = encode.se, labels = encode.se$label.main)
s2022_so$SingleR_hpca <- predicted_cell_types_HPCA$labels
s2022_so$SingleR_encode <- predicted_cell_types_ENCODE$labels


## plot metadata
to_plot <- c('sample', 'patient', 'tissue', 'timepoint', 'leukemia', 
             'subtype', 'relapse', 'Phase', 'SingleR_hpca', 'SingleR_encode')

for(val in to_plot){
  name_nolab <- paste0('s2022_images/DimPlot_', val, '_nolab.pdf')
  cairo_pdf(name_nolab, width = 6, height = 6)
  print(DimPlot(s2022_so, group.by = val) & NoLegend())
  dev.off()
  name_lab <- paste0('s2022_images/DimPlot_', val, '.pdf')
  cairo_pdf(name_lab, width = 8, height = 6)
  print(DimPlot(s2022_so, group.by = val))
  dev.off()
}

## save Seurat object
saveRDS(s2022_so, 'SeuratObjects/s2022_processed.rds')


### Merge Seurat objects ----- 

DefaultAssay(pilot_so) <- "SCT"
DefaultAssay(s2022_so) <- "SCT"

so <- merge(pilot_so, s2022_so, merge.data = TRUE)
so
#An object of class Seurat 
#69230 features across 144747 samples within 2 assays 
#Active assay: SCT (32629 features, 0 variable features)
#3 layers present: counts, data, scale.data
#1 other assay present: RNA

#set variable features
#per https://github.com/satijalab/seurat/issues/2814
VariableFeatures(so[["SCT"]]) <- rownames(so[["SCT"]]@scale.data)

## Dimensionality reduction & clustering
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)

so <- FindNeighbors(so, dims = 1:30)
so <- FindClusters(so) 

## plot metadata
to_plot <- c('batch', 'sample', 'patient', 'tissue', 'timepoint', 'leukemia', 
             'subtype', 'relapse', 'Phase', 'SingleR_hpca', 'SingleR_encode')

for(val in to_plot){
  name_nolab <- paste0('merge_images/DimPlot_', val, '_nolab.pdf')
  cairo_pdf(name_nolab, width = 6, height = 6)
  print(DimPlot(so, group.by = val) & NoLegend())
  dev.off()
  name_lab <- paste0('merge_images/DimPlot_', val, '.pdf')
  cairo_pdf(name_lab, width = 8, height = 6)
  print(DimPlot(so, group.by = val))
  dev.off()
}

## save Seurat object
saveRDS(so, 'SeuratObjects/merge_all.rds')

## plot composition

so$timepoint_leukemia <- paste(so$timepoint, so$leukemia, sep = '_')

groups <- c('sample', 'patient', 'leukemia', 'timepoint_leukemia')
compositions <- c('Phase', 'SingleR_hpca', 'SingleR_encode')

for(group in groups){
  for(composition in compositions){
    data <- FetchData(so, vars = c(group, composition))
    table <- table(data[[1]], data[[2]])
    name <- paste0('merge_images/composition_', group, '_', composition, '.pdf')
    width <- case_when(
      group == 'sample' ~ 12,
      group == 'patient' ~ 10,
      TRUE ~ 6
    )
    cairo_pdf(name, width = width, height = 6)
    print(ggplot(as.data.frame(t(table)), aes(fill = Var1, y = Freq, x = Var2)) + 
            geom_bar(position = 'fill', stat = 'identity') + RotatedAxis())
    dev.off()
  }
}


## leukemia DimPlot redo to use leukemia color scheme

leukemia_color_list <- hue_pal()(7)
names(leukemia_color_list) <- unique(sampleTable$leukemia)

cairo_pdf('merge_images/DimPlot_leukemia_nolab.pdf', width = 6, height = 6)
print(DimPlot(so, group.by = 'leukemia', cols = leukemia_color_list) & NoLegend())
dev.off()
cairo_pdf('merge_images/DimPlot_leukemia.pdf', width = 8, height = 6)
print(DimPlot(so, group.by = 'leukemia', cols = leukemia_color_list))
dev.off()


