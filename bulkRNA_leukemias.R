library(tidyverse)
library(pheatmap)
library(DESeq2)
library(limma)
library(scales)
library(factoextra)

sampleTable <- read.csv('sampleTable.csv')
rownames(sampleTable) <- sampleTable$ID

## define color scheme for leukemias
length(unique(sampleTable$leukemia))
# 7

leukemia_color_list <- hue_pal()(7)
names(leukemia_color_list) <- unique(sampleTable$leukemia)

#create color scale
leukemia_col_scale <- scale_color_manual(name = "leukemia", values = leukemia_color_list)


#filter down to bulk RNAseq samples only
sampleTable_bulk <- sampleTable[sampleTable$b_batch != '', ]

raw_counts <- read.csv('inputs/bulkRNA_kallisto_rawCounts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(raw_counts) <- gsub('\\.', '-', colnames(raw_counts))


### Use DESeq2 and limma to normalize and batch correct -----

#set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(raw_counts), 
                              colData = sampleTable_bulk,
                              design = ~ b_batch + leukemia)

#calculate size factors for normalization
dds <- estimateSizeFactors(dds)
#log2(normalized counts + 1) with pseudocount
log2_counts <- log2(counts(dds, normalized = TRUE) + 1)

## batch correct using limma
design <- model.matrix(~leukemia, sampleTable_bulk)
log2_counts_bc <- limma::removeBatchEffect(log2_counts, sampleTable_bulk$b_batch, design = design)

#save these
write.csv(log2_counts_bc, 'all_bulk_log2normalized_bc_counts.csv', row.names = TRUE, quote = FALSE)


### Filter down to diagnosis (DX) samples only -----

sampleTable_DX <- sampleTable_bulk[sampleTable_bulk$timepoint == 'Diagnosis', ]

table(sampleTable_DX$leukemia)
#AMKL   AML B-ALL   MDS T-ALL 
#   1     9    31     1     7 

counts_DX <- raw_counts[ , colnames(raw_counts) %in% sampleTable_DX$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(counts_DX), 
                              colData = sampleTable_DX,
                              design = ~ b_batch + leukemia)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~leukemia, sampleTable_DX)
mat <- limma::removeBatchEffect(mat, vsd$b_batch, design = design0)
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
pdf('bulkRNA_leukemias/PCA.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_DX$leukemia, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 5)) + leukemia_col_scale
dev.off()
# PCA plot colored according to tissue type (same patient samples connected)
pdf('bulkRNA_leukemias/PCA_tissue.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_DX$tissue, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 5)) + 
  geom_line(aes(group = sampleTable_DX$patient), color = 'grey30', linewidth = 2)
dev.off()


### Differential expression -----

dds <- DESeq(dds, test = "LRT", reduced = ~ b_batch)

res <- results(dds, alpha = 0.05)

summary(res)
#out of 23280 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 3848, 17%
#LFC < 0 (down)     : 2772, 12%
#outliers [1]       : 247, 1.1%
#low counts [2]     : 2687, 12%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#convert to dataframe
res_df <- as.data.frame(res)

#filter down to significantly DE genes
sig_res <- res_df[res_df$padj < 0.05, ]
sig_res <- sig_res[!is.na(sig_res$padj), ]
sig_res <- sig_res[order(sig_res$pad), ]

write.csv(sig_res, 'bulkRNA_leukemias/DE_sig_genes.csv', row.names = TRUE, quote = FALSE)


### Plot top DE genes -----

#read in normalized counts
log_counts <- read.csv('all_bulk_log2normalized_bc_counts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(log_counts) <- gsub('\\.', '-', colnames(log_counts))

#filter down to DX samples only
log_counts_DX <- log_counts[ , colnames(log_counts) %in% sampleTable_DX$ID]

#pull out genes of interest
goi <- rownames(sig_res[1:50, ])

#filter counts to goi
coi <- log_counts_DX[rownames(log_counts_DX) %in% goi, ]

color_list <- list(leukemia = leukemia_color_list[names(leukemia_color_list) %in% unique(sampleTable_DX$leukemia)])

pdf('bulkRNA_leukemias/top50DE_heatmap.pdf', height = 8, width = 12)
pheatmap(t(coi), scaling = 'none', 
         annotation_row = sampleTable_DX['leukemia'], 
         annotation_colors = color_list)
dev.off()


### Compare to FDA Pediatric Molecular Target Lists (https://moleculartargets.ccdi.cancer.gov/fda-pmtl) -----

PMTL <- read.csv('FDA_PediatricMolecularTargets.csv')

#find overlap with differentially expressed genes
goi_PMTL <- sig_res[rownames(sig_res) %in% PMTL$targetSymbol, ]
#160 genes

#plot top 10
coi_PMTL <- log_counts_DX[rownames(log_counts_DX) %in% rownames(goi_PMTL)[1:10], ]
t_coi_PMTL <- as.data.frame(t(coi_PMTL))

#create fill scale
leukemia_fill_scale <- scale_fill_manual(name = "leukemia", values = leukemia_color_list)

for(gene in colnames(t_coi_PMTL)){
  filename <- paste0('bulkRNA_leukemias/geneexp_', gene, '.pdf')
  pdf(filename, height = 4, width = 3)
  print(ggplot(t_coi_PMTL, aes(x = sampleTable_DX$leukemia, y = .data[[gene]], 
                               fill = sampleTable_DX$leukemia)) + 
          geom_boxplot() + geom_jitter(width = 0.1, height = 0.1) + 
          leukemia_fill_scale + xlab(gene) + ylab('log(counts)') + 
          theme_classic() & RotatedAxis())
  dev.off()
}





