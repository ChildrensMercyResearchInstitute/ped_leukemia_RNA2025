library(tidyverse)
library(pheatmap)
library(DESeq2)
library(limma)
library(scales)
library(factoextra)
library(gprofiler2)
library(Polychrome)
library(ggvenn)

sampleTable <- read.csv('inputs/sampleTable.csv')
rownames(sampleTable) <- sampleTable$ID

## using color scheme for subtypes defined in Fig 1
subtype_color_list <- kelly.colors(20)[2:20]
names(subtype_color_list) <- unique(sampleTable$subtype)

#create color scale
subtype_col_scale <- scale_color_manual(name = "subtype", values = subtype_color_list)

## filter down to B-ALL samples at diagnosis only

sampleTable_BALL <- sampleTable[sampleTable$leukemia == 'B-ALL', ]
sampleTable_BALL <- sampleTable_BALL[sampleTable_BALL$timepoint == 'Diagnosis', ]


#separate bulk and sc samples
sampleTable_bulk <- sampleTable_BALL[sampleTable_BALL$b_batch != '', ]
sampleTable_sc <- sampleTable_BALL[sampleTable_BALL$sc_batch != '', ]

table(sampleTable_bulk$subtype)
#   BCR::ABL1-like    ETV6::RUNX1   hyperdiploid      IGH::BCL2        KMT2A-r       ZNF384-r 
# 6              6             11              5              1              1              1 

table(sampleTable_sc$subtype)
#        BCR::ABL1 BCR::ABL1-like    ETV6::RUNX1   hyperdiploid    hypodiploid        KMT2A-r     TCF3::PBX1 
# 3              1              7             14              5              1              2              1 
# ZNF384-r 
#        1 


### bulk visualization -----
#read in counts
raw_counts <- read.csv('inputs/bulkRNA_kallisto_rawCounts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(raw_counts) <- gsub('\\.', '-', colnames(raw_counts))

bulk_counts <- raw_counts[ , colnames(raw_counts) %in% sampleTable_bulk$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(bulk_counts), 
                              colData = sampleTable_bulk,
                              design = ~ b_batch + subtype)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~subtype, sampleTable_bulk)
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

# PCA plot colored according to subtype
pdf('B-ALL_subtypes/bulk_PCA.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_bulk$subtype, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 7)) + subtype_col_scale
dev.off()


### single-cell pseudobulk visualization -----
#read in counts
pbulk_counts <- read.csv('inputs/pseudobulk_sampleCounts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(pbulk_counts) <- gsub('\\.', '-', colnames(pbulk_counts))

sc_counts <- pbulk_counts[ , colnames(pbulk_counts) %in% sampleTable_sc$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(sc_counts), 
                              colData = sampleTable_sc,
                              design = ~ sc_batch + subtype)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~subtype, sampleTable_sc)
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

# PCA plot colored according to subtype
pdf('B-ALL_subtypes/sc_PCA.pdf', width = 8, height = 6)
fviz_pca_ind(pca, habillage = sampleTable_sc$subtype, geom = "point", 
             mean.point = FALSE, pointsize = 4) + 
  scale_shape_manual(values = rep(16, 9)) + subtype_col_scale
dev.off()


### bulk differential expression across known subtypes -----

sampleTable_bulk_known <- sampleTable_bulk[sampleTable_bulk$subtype != '', ]
bulk_counts_known <- bulk_counts[ , colnames(bulk_counts) %in% sampleTable_bulk_known$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(bulk_counts_known), 
                              colData = sampleTable_bulk_known,
                              design = ~ b_batch + subtype)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~subtype, sampleTable_bulk_known)
mat <- limma::removeBatchEffect(mat, vsd$b_batch, design = design0)
assay(vsd) <- mat

## differential expression
dds <- DESeq(dds, test = "LRT", reduced = ~ b_batch)

res <- results(dds, alpha = 0.05)

summary(res)
#out of 22880 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2133, 9.3%
#LFC < 0 (down)     : 3442, 15%
#outliers [1]       : 176, 0.77%
#low counts [2]     : 3519, 15%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#convert to dataframe
bulk_df <- as.data.frame(res)

#filter down to significantly DE genes
sig_bulk <- bulk_df[bulk_df$padj < 0.05, ]
sig_bulk <- sig_bulk[!is.na(sig_bulk$padj), ]
sig_bulk <- sig_bulk[order(sig_bulk$pad), ]

write.csv(sig_bulk, 'B-ALL_subtypes/bulk_DE_sig_genes.csv', row.names = TRUE, quote = FALSE)


### sc differential expression across known subtypes -----

sampleTable_sc_known <- sampleTable_sc[sampleTable_sc$subtype != '', ]
#due to transcriptional similarity between BCR::ABL1 and BCR::ABL1-like, 
#consider the single BCR::ABL1 case as BCR::ABL1-like for DE purposes
sampleTable_sc_known$subtype <- case_when(
  sampleTable_sc_known$subtype == 'BCR::ABL1' ~ 'BCR::ABL1-like',
  TRUE ~ sampleTable_sc_known$subtype
)

sc_counts_known <- sc_counts[ , colnames(sc_counts) %in% sampleTable_sc_known$ID]

## set up DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(sc_counts_known), 
                              colData = sampleTable_sc_known,
                              design = ~ sc_batch + subtype)

#transform counts for data visualization
vsd <- vst(dds, blind = FALSE)
#batch correction
mat <- assay(vsd)
design0 <- model.matrix(~subtype, sampleTable_sc_known)
mat <- limma::removeBatchEffect(mat, vsd$sc_batch, design = design0)
assay(vsd) <- mat

## differential expression
dds <- DESeq(dds, test = "LRT", reduced = ~ sc_batch)

res <- results(dds, alpha = 0.05)

summary(res)
#out of 34114 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1089, 3.2%
#LFC < 0 (down)     : 1973, 5.8%
#outliers [1]       : 133, 0.39%
#low counts [2]     : 16407, 48%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#convert to dataframe
sc_df <- as.data.frame(res)

#filter down to significantly DE genes
sig_sc <- sc_df[sc_df$padj < 0.05, ]
sig_sc <- sig_sc[!is.na(sig_sc$padj), ]
sig_sc <- sig_sc[order(sig_sc$pad), ]

write.csv(sig_sc, 'B-ALL_subtypes/sc_DE_sig_genes.csv', row.names = TRUE, quote = FALSE)


### explore overlap (including with PMTL) -----
overlap <- sig_bulk[rownames(sig_bulk) %in% rownames(sig_sc), ]
#1216 genes

PMTL <- read.csv('FDA_PediatricMolecularTargets.csv')

overlap_all <- overlap[rownames(overlap) %in% PMTL$targetSymbol, ]
#41 genes

write.csv(overlap_all, 'B-ALL_subtypes/subtype_bulk_sc_PMTL_overlapGenes.csv', row.names = TRUE, quote = FALSE)

## Venn diagram
venn_sets <- list(bulkRNA = unique(rownames(sig_bulk)), 
                  scRNA = unique(rownames(sig_sc)), 
                  PMTL = unique(PMTL$targetSymbol))
pdf('B-ALL_subtypes/DE_PMTL_Venn.pdf', height = 4, width = 4)
ggvenn(venn_sets, show_percentage = FALSE, 
       fill_color = c('white', 'white', 'white'))
dev.off()

### visualize expression across bulk RNAseq dataset -----

#read in normalized, batch-corrected counts
log_counts <- read.csv('all_bulk_log2normalized_bc_counts.csv', row.names = 1)
#replace '.' in colnames with '-' to match sampleTable
colnames(log_counts) <- gsub('\\.', '-', colnames(log_counts))

#filter counts to only B-ALL samples
log_counts_BALL <- log_counts[ , colnames(log_counts) %in% sampleTable_bulk$ID]

#pull out genes of interest (top 500)
goi <- rownames(overlap[1:500, ])

#filter counts to genes of interest
coi <- log_counts_BALL[rownames(log_counts_BALL) %in% goi, ]

color_list <- list(subtype = subtype_color_list[names(subtype_color_list) %in% unique(sampleTable_bulk$subtype)])

pdf('B-ALL_subtypes/bulk_DEgenes_heatmap.pdf', height = 8, width = 12)
pheatmap(t(coi), scaling = 'none', show_colnames = FALSE,
         annotation_row = sampleTable_bulk['subtype'], 
         annotation_colors = color_list)
dev.off()



