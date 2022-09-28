library(readxl)
library(writexl)
library(tidyverse)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)


cData <- readxl::read_xlsx(
    path = './../data/macrogen_mapping_result.xlsx',
    col_names = TRUE
)

cData <- cData %>% column_to_rownames(var = "Gene_Symbol")
colnames(cData) <- sapply(colnames(cData), function(x) str_split(x, '_')[[1]][1])

# margin == 1 : row-wise
cData <- cData[
    apply(cData, 1, function(x) all(x != 0)), 
]

# save original count data
raw_cData <- cData

cData <- cData[, -c(4, 5)]


mData <- read.table(
    file = './../data/meta_data',
    header = TRUE
)

# save original meta data
raw_mData <- mData

raw_mData <- raw_mData %>% column_to_rownames(var = "SampleName")

raw_mData$SampleType <- as.factor(raw_mData$SampleType)
raw_mData$SampleType <- relevel(raw_mData$SampleType, ref = 'CON')


mData <- mData %>% as.data.frame()
mData <- mData[-c(4, 5), ]
rownames(mData) <- NULL

mData <- mData %>% column_to_rownames(var = "SampleName")

mData$SampleType <- as.factor(mData$SampleType)
mData$SampleType <- relevel(mData$SampleType, ref = 'CON')


dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cData,
    colData = mData,
    design = ~SampleType
)

ddsX <- DESeq2::DESeq(dds)


rld_T <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld_T)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], SampleType = mData)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = SampleType, size = 5)) +
    xlab(paste0("PC1 : ", percentVar[1], "%")) +
    ylab(paste0("PC2 : ", percentVar[2], "%")) +
    coord_fixed(ratio = sd_ratio) +
    geom_text_repel(aes(label = rownames(pca_df)))


rld_cor <- stats::cor(rld_mat)

heat_colors <- rev(RColorBrewer::brewer.pal(11, 'RdBu'))
pheatmap::pheatmap(mat = rld_cor, annotation = mData, color = heat_colors)


lfcR <- DESeq2::lfcShrink(ddsX, 
                          coef = DESeq2::resultsNames(ddsX)[2], 
                          type = 'apeglm')

sig_lfcR <- lfcR %>% 
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "Gene")

sig_lfcR2 <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    rownames_to_column(var = "Gene")

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(Gene %in% sig_lfcR$Gene)

sig_norm_cnt2 <- norm_cnt %>%
    dplyr::filter(Gene %in% sig_lfcR2$Gene)


pheatmap::pheatmap(
    mat = sig_norm_cnt %>% column_to_rownames(var = "Gene"),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = 'row',
    color = heat_colors
)

pheatmap::pheatmap(
    mat = sig_norm_cnt2 %>% column_to_rownames(var = "Gene"),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = 'row',
    color = heat_colors
)


# save DEG result
writexl::write_xlsx(
    x = sig_lfcR,
    path = './../results/sig_lfcR.xlsx',
    col_names = TRUE
)

writexl::write_xlsx(
    x = sig_lfcR2,
    path = './../results/sig_lfcR2.xlsx',
    col_names = TRUE
)


# Plot top 20 Down-regulated DEGs
down_sig <- readxl::read_excel(
    path = './../results/sig_lfcR_Down.xlsx',
    col_names = TRUE
)

down_sig <- down_sig[order(down_sig$log2FoldChange, decreasing = FALSE), ]

## Gene names should not be ordered alphabetically --> order with log2FoldChange
## in here, we already ordered gene name by log2FoldChange above
## so, we just plugin (already ordered) gene names to levels parameter
down_sig$Gene <- factor(down_sig$Gene, levels = rev(down_sig$Gene))

ggplot(down_sig[1:30, ], aes(x = -log2FoldChange, y = Gene)) + 
    geom_point(aes(size = 1, color = -log(padj))) +
    scale_color_gradient(low = "#1E88E5", high = "#E53935") +
    xlim(1, 4.5)

ggplot(down_sig[1:30, ], 
       aes(x = -log2FoldChange, y = Gene, fill = -log(padj))) + 
    geom_bar(stat = 'identity') +
    scale_fill_gradient(low = "#1E88E5", high = "#E53935")


# GO and KEGG result
func_res <- readxl::read_excel(
    path = './../results/DAVID_GO_KEGG_renamed.xlsx',
    col_names = TRUE
)

func_res <- func_res[order(func_res$Count, decreasing = FALSE), ]

func_res$Term <- factor(func_res$Term, levels = func_res$Term)

# Plot KEGG & GO analysis result
#@@ 1. Bar plot & Dot plot
ggplot(func_res, aes(x = Count, y = Term)) + 
    geom_point(aes(size = 1, color = -log(FDR))) +
    scale_color_gradient(low = "#1E88E5", high = "#E53935") +
    xlim(0, 40)

ggplot(func_res, aes(x = Count, y = Term, fill = -log(FDR))) + 
    geom_bar(stat = 'identity', width = 0.3) +
    scale_fill_gradient(low = "#1E88E5", high = "#E53935")

### 2. 