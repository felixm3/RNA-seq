# Suppress warnings
options(warn = -1)

# load packages
library(tximport)
library(ensembldb)
library(AnnotationHub)
library(DESeq2)

library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

library(biomaRt)

library(EnhancedVolcano)

# increase size of plots from default
options(repr.plot.width = 14, 
        repr.plot.height = 14) # from 7, 7

# create an AnnotationHub object
ah = AnnotationHub()

# find relevant annotation records
ahs2 <- query(ah, c("110","Ensembl", "Mus musculus"))

# load latest annotation record for Mus musculus
ensdb_110 <- ahs2[["AH113713"]]

# Extract transcript and gene information
tx_data <- transcripts(ensdb_110, return.type = "DataFrame")

# Create the tx2gene data.frame
tx2gene <- tx_data[, c("tx_id", "gene_id")]
tx2gene

# prep salmon quant files for the 4 samples for tximport
quants_dir <- "results/03_salmon_quant/"
all_files <- list.files(quants_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
quant_files <- all_files[grepl("SRR.*quant.sf", all_files)] # Use grepl to select only quant.sf files
quant_dirs <- list.files(quants_dir, pattern = "quant_SRR.*", full.names = TRUE)
sample_names <- sub('.*SRR([0-9]+).*', 'SRR\\1', quant_dirs)  # get sample names from dir names
names(quant_files) <- sample_names
as.data.frame(quant_files)

# tximport quant_files
txi <- tximport(quant_files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)
str(txi)

# create sample table
condition <- factor(c('young', 'young', 'aged', 'aged'))
coldata <- data.frame(row.names = sample_names, condition)
coldata

# create deseq object
dds <- DESeqDataSetFromTximport(txi, coldata, ~ condition)
dds

# drop genes with fewer than 12 counts total across all samples
keep <- rowSums(counts(dds)) > 12
dds <- dds[keep, ]
dds

# rlog transformation to stabilize variance
rld <- rlog(dds)

# sample distances
sampleDists <- dist(t(assay(rld)))

# heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, colnames(rld), sep = " - " )
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, 
        angle_col = 45, fontsize = 20, 
        display_numbers = TRUE)

# PCA
object <- rld
ntop <- 500 # number of variable genes to use for PCA
intgroup <- 'condition'
pcsToUse = 1:2

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

pcs <- paste0("PC", pcsToUse)
d <- data.frame(V1 = pca$x[, pcsToUse[1]], V2 = pca$x[, pcsToUse[2]], 
                group = group, intgroup.df, name = colnames(object))
colnames(d)[1:2] <- pcs

ggplot(data = d, aes(x = PC1, y = PC2, color = group)) + 
    geom_point(size = 16) + 
    xlab(paste0(pcs[1], ": ", round(percentVar[pcsToUse[1]] * 100), "% variance")) + 
    ylab(paste0(pcs[2], ": ", round(percentVar[pcsToUse[2]] * 100), "% variance")) + 
    coord_fixed() + 
    theme_gray(base_size = 24)

dds <- DESeq(dds)
dds

res <- results(dds)
summary(res)
res

# using more stringent thresholds: BH FDR adjusted p-value 0.01 vs 0.1, LFC 0.1 vs 0
res <- results(dds, alpha = 0.01, lfcThreshold=0.1)
summary(res)
res

# add BIOMART gene symbols to res based on Ensembl IDs in rownames
ensembl <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
res$symbol <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                          filters = "ensembl_gene_id",
                          values = rownames(res),
                          mart = ensembl)$external_gene_name
res

# what genes have the strongest down-regulation in young?
resSig <- subset(res, padj < 0.01)
head(resSig[ order(resSig$log2FoldChange), ])

# what genes have the strongest UP-regulation in young?
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

EnhancedVolcano(as.data.frame(res), x = 'log2FoldChange', y = 'padj', lab = as.data.frame(res)$'symbol', 
               pCutoff = 1e-06, FCcutoff = 3, 
               xlim = c(-12, 12))

# heatmap of genes with strongest signal
resSig <- subset(res, padj < 1e-6 & abs(log2FoldChange) > 3)
# Create a named vector for mapping Ensembl IDs to gene symbols
ensembl_to_symbol <- setNames(res$symbol, rownames(res))
# select only resSig genes instead of heatmap of all ~17,000 genes
mat  <- assay(rld)[rownames(resSig), ] 
# heatmap of deviation from mean expression instead of actual expression
mat  <- mat - rowMeans(mat) 
# annotate heatmap columns with sample name
anno <- as.data.frame(colData(rld))
pheatmap(mat, 
         annotation_col = anno, 
        labels_row = ensembl_to_symbol[rownames(mat)], 
        fontsize = 20)

resSig[ order(resSig$padj, decreasing = FALSE), c('baseMean', 'log2FoldChange', 'padj', 'symbol')]

sessionInfo()

