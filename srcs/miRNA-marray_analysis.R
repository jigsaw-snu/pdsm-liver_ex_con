# [Acknowledgement]
# - Affymetrix array : single channel (one color)
# - Agilent array : Single / Dual channel (One or Two color)
# - Illumina array : Single channel (one color)
#
#
# [References]
# - Microarray Analysis Workflow (from CEL to Differential Analysis)
#   http://homer.ucsd.edu/homer/basicTutorial/affymetrix.html
#   https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
#
# - Explains for RMA function in oligo package
#   https://gtk-teaching.github.io/Microarrays-R/05-DataNormalisation/index.html


# install necessary pacakges
# Try use Biocmanager::install()
# If installation not available, then use install.packages() to install from source
install.packages("~/Downloads/oligo_1.58.0.tar.gz", 
                 repos = NULL, 
                 type = "source")

install.packages("~/Downloads/pd.mirna.4.0_3.12.0.tar.gz",
                 repos = NULL,
                 type = "source")


library(Biobase)  # General Bioconductor package ; i.e. pData(), fData(), exprs()
#library(affy) # read CEL files for old version of Affymetrix arrays
library(oligo)  # read CEL files ; only for Affymetrix and Nimblegen microarrys
library(pd.mirna.4.0)  # probe annotation (package depends on type of array platform)
library(biomaRt)
#library(mirbase.db)
library(miRBaseConverter)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(plotly)  # 3d scatter plot --> plot PCA in 3-Dimensional space


# reading raw data files (CEL files)
raw_file_names <- list.files('./../data/CEL_Files/', full.names = TRUE)

# raw_file_names <- raw_file_names[-c(2, 6)]

raw_data <- oligo::read.celfiles(
    filenames = raw_file_names
)

# modify sample names and add extra information to pData
colnames(raw_data) <- c(
    'CON_Liver-1', 'CON_Liver-2', 'CON_Liver-5', 'EX_Liver-2',
    'EX_Liver-3', 'EX_Liver-5', 'EX_Liver-6', 'CON_Liver-4'
)

# colnames(raw_data) <- c(
#     'CON_Liver-1', 'CON_Liver-5', 'EX_Liver-2',
#     'EX_Liver-3', 'EX_Liver-6', 'CON_Liver-4'
# )

Biobase::pData(raw_data)$Sample.Group <- c(rep('CON', 3), rep('EX', 4), 'CON')

# Biobase::pData(raw_data)$Sample.Group <- c(rep('CON', 2), rep('EX', 3), 'CON')

# RMA rules over any other normalization algorithms in microarray analysis thesedays
# RMA algorithm performs
# 1. Background correction : remove backgound noise
# 2. Normalization : quantile normalization of probe level data
# 3. Summarization : median polish algorithm -> compress expr data into a gene/probeset data
eset <- oligo::rma(raw_data)


# QC : plot PCA
pca <- stats::prcomp(t(Biobase::exprs(eset)), scale. = TRUE)  # transpose is needed to get PC of each 'sample', not gene(probe)

percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)  # it seems top loadings are similar

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
    PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6],
    Group = Biobase::pData(eset)$Sample.Group
)

# pca-2D
ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = Group, size = 5)) + 
    xlab(paste0('PC1 ', 'VarExp: ', percentVar[1], '%')) + 
    ylab(paste0('PC2 ', 'VarExp: ', percentVar[2], '%')) +
    coord_fixed(ratio = sd_ratio) +
    geom_text_repel(aes(label = rownames(pca_df)))


# scree plot
qplot(x = c(1:8), y = percentVar) + 
    geom_line() +
    xlab('Principal Component') +
    ylab('Variance Explained') +
    ylim(0, 16)


# pca-3D
plotly::plot_ly(pca_df,
                x = pca_df$PC1, y = pca_df$PC2, z = pca_df$PC3,
                color = pca_df$Group, colors = c('#5DADE2', '#E74C3C'))



# QC : plot Heatmap with cluster dendrogram
cor <- stats::cor(Biobase::exprs(eset))

heatcolor <- rev(RColorBrewer::brewer.pal(11, 'PRGn'))

pheatmap::pheatmap(cor, color = heatcolor)




# Differential Expression Analysis
library(limma)  # differential expression analysis

groups <- Biobase::pData(eset)$Sample.Group
des_fct <- factor(groups, levels = c('CON', 'EX'))  # base level = CON
limma_design <- stats::model.matrix(~ 0 + des_fct)
colnames(limma_design) <- c('CON', 'EX')

fit <- limma::lmFit(Biobase::exprs(eset), limma_design)

cont <- limma::makeContrasts(EX-CON, levels = limma_design)
fit.cont <- limma::contrasts.fit(fit, cont)

fit.eb <- limma::eBayes(fit.cont)

limma_res <- limma::topTable(fit.eb, number = Inf, lfc = 0.58, adjust.method = 'BH')

sig_limma <- limma_res %>% dplyr::filter(adj.P.Val < 0.05)

test_sig_lm <- limma::topTable(fit.eb, 
                             number = Inf,
                             lfc = 1, 
                             adjust.method = 'BH') %>% 
    dplyr::filter(P.Value < 0.05)

limma_res_with_annot <- limma::topTable(fit.eb, 
                                        number = Inf, 
                                        lfc = 0, 
                                        adjust.method = 'BH')
limma_res_with_annot$reg <- rep('None', nrow(limma_res_with_annot))
limma_res_with_annot$id <- rep('', nrow(limma_res_with_annot))

for (i in 1:nrow(limma_res_with_annot)) {
    if (limma_res_with_annot$adj.P.Val[i] < 0.05) {
        if (limma_res_with_annot$logFC[i] < -0.58) {
            limma_res_with_annot$reg[i] <- 'down'
            limma_res_with_annot$id[i] <- rownames(limma_res_with_annot)[i]
        }
        
        if (limma_res_with_annot$logFC[i] > 0.58) {
            limma_res_with_annot$reg[i] <- 'up'
            limma_res_with_annot$id[i] <- rownames(limma_res_with_annot)[i]
        }
    }
}

# volcano plot
ggplot(limma_res_with_annot,
       aes(x = logFC, y = -log(adj.P.Val), color = reg)) + 
    geom_point() +
    scale_color_manual(
        values = c('down' = 'blue', 'up' = 'red', 'None' = 'grey')
    ) + 
    geom_text_repel(aes(label = id)) +
    ylim(0, 4)


# Get Gene Symbols first
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
dset <- biomaRt::useDataset(mart = mart,
                            dataset = 'mmusculus_gene_ensembl')