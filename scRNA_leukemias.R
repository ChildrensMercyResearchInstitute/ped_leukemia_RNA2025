library(tidyverse)
library(pheatmap)
library(DESeq2)
library(limma)
library(scales)
library(factoextra)
library(gprofiler2)

sampleTable <- read.csv('inputs//sampleTable.csv')
rownames(sampleTable) <- sampleTable$ID

## define color scheme for leukemias
length(unique(sampleTable$leukemia))
# 7

leukemia_color_list <- hue_pal()(7)
names(leukemia_color_list) <- unique(sampleTable$leukemia)

#create color scale
leukemia_col_scale <- scale_color_manual(name = "leukemia", values = leukemia_color_list)


#filter down to single-cell RNAseq samples only
sampleTable_sc <- sampleTable[sampleTable$sc_batch != '', ]

#read in pseudobulk counts
pbulk_counts <- read.csv('inputs/pseudobulk_sampleCounts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(pbulk_counts) <- gsub('\\.', '-', colnames(pbulk_counts))


### Use DESeq2 and limma to normalize and batch correct -----

#set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(pbulk_counts), 
                              colData = sampleTable_sc,
                              design = ~ sc_batch + leukemia)

#calculate size factors for normalization
dds <- estimateSizeFactors(dds)
#log2(normalized counts + 1) with pseudocount
log2_counts <- log2(counts(dds, normalized = TRUE) + 1)

## batch correct using limma
design <- model.matrix(~leukemia, sampleTable_sc)
log2_counts_bc <- limma::removeBatchEffect(log2_counts, sampleTable_sc$sc_batch, design = design)

#save these
write.csv(log2_counts_bc, 'all_sc_pseudobulk_log2normalized_bc_counts.csv', row.names = TRUE, quote = FALSE)


### Filter down to diagnosis (DX) samples only -----

sampleTable_DX <- sampleTable_sc[sampleTable_sc$timepoint == 'Diagnosis', ]

table(sampleTable_DX$leukemia)
#ALCL  AMKL   AML B-ALL   MDS  MPAL T-ALL 
#   1     1    13    35     1     2    11 

counts_DX <- pbulk_counts[ , colnames(pbulk_counts) %in% sampleTable_DX$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(counts_DX), 
                              colData = sampleTable_DX,
                              design = ~ sc_batch + leukemia)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~leukemia, sampleTable_DX)
mat <- limma::removeBatchEffect(mat, vsd$sc_batch, design = design0)
assay(vsd) <- mat

## manually run PCA for improved plotting
# code adapted from DESeq2 plotPCA function

# calculate the variance for each gene
rv <- rowVars(assay(vsd))
# select the 500 genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# PCA plot colored according to leukemia
pdf('scRNA_leukemias/PCA.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_DX$leukemia, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 7)) + leukemia_col_scale
dev.off()
# PCA plot colored according to tissue type (same patient samples connected)
pdf('scRNA_leukemias/PCA_tissue.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_DX$tissue, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 7)) + 
  geom_line(aes(group = sampleTable_DX$patient), color = 'grey30', linewidth = 2)
dev.off()


### Differential expression -----

dds <- DESeq(dds, test = "LRT", reduced = ~ sc_batch)

res <- results(dds, alpha = 0.05)

summary(res)
#out of 34770 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2839, 8.2%
#LFC < 0 (down)     : 1956, 5.6%
#outliers [1]       : 148, 0.43%
#low counts [2]     : 14084, 41%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#convert to dataframe
res_df <- as.data.frame(res)

#filter down to significantly DE genes
sig_res <- res_df[res_df$padj < 0.05, ]
sig_res <- sig_res[!is.na(sig_res$padj), ]
sig_res <- sig_res[order(sig_res$pad), ]

write.csv(sig_res, 'scRNA_leukemias/DE_sig_genes.csv', row.names = TRUE, quote = FALSE)


### Look for overlap with bulk RNAseq DE -----

bulk_DE <- read.csv('bulkRNA_leukemias/DE_sig_genes.csv', row.names = 1)

overlap <- sig_res[rownames(sig_res) %in% rownames(bulk_DE), ]
#1958 genes

## pathway analysis

gostres <- gost(query = rownames(overlap), organism = "hsapiens", ordered_query = FALSE, 
                sources = c("GO:BP", "GO:MF", "KEGG", "REAC", "CORUM"),
                multi_query = FALSE, significant = TRUE, user_threshold = 0.05, 
                correction_method = "g_SCS", domain_scope = "annotated")

pathways_df <- as.data.frame(gostres$result)

#save (without 'parents' column which is a list)
write.csv(pathways_df[1:13], 'scRNA_leukemias/bulk_sc_overlap_pathways.csv', row.names = FALSE, quote = TRUE)

## overlap with PMTL

PMTL <- read.csv('FDA_PediatricMolecularTargets.csv')

overlap_all <- overlap[rownames(overlap) %in% PMTL$targetSymbol, ]
#66 genes

write.csv(overlap_all, 'scRNA_leukemias/bulk_sc_PMTL_overlapGenes.csv', row.names = TRUE, quote = FALSE)
