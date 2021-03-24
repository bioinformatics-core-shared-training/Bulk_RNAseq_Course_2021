library(EnsDb.Mmusculus.v79)
library(DESeq2)
library(tidyverse)

# load the deseq2 results

ddsObj.interaction <- readRDS("RObjects/DESeqDataSet.interaction.rds")
results.interaction.11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
results.interaction.33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")

results.interaction.11

# Annotation

columns(EnsDb.Mmusculus.v79)

keytypes(EnsDb.Mmusculus.v79)

ourCols <- c("SYMBOL", "GENEID", "ENTREZID")
ourKeys <- rownames(results.interaction.11)[1:1000]

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
annot
dim(annot)
length(unique(annot$ENTREZID))

# Exercise 1 

ourKeys <- rownames(results.interaction.11)
ourCols <- c("SYMBOL", "GENEID", "ENTREZID", "GENEBIOTYPE")

annot <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                               keys = ourKeys,
                               columns = ourCols,
                               keytype = "GENEID")
length(ourKeys)
length(unique(annot$GENEID))

# One we made earlier

ensemblAnnot <- readRDS("RObjects/Ensembl_annotations.rds")

annot.interaction.11 <- as.data.frame(results.interaction.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC = log2FoldChange, FDR = padj)

# Visualisation

hist(annot.interaction.11$pvalue)

# shrinking the log 2 fold change

ddsShrink.11 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.11,
                          type = "ashr")
shrinkTab.11 <- as.data.frame(ddsShrink.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC = log2FoldChange, FDR = padj)

# MA plots

plotMA(results.interaction.11, alpha = 0.05)
plotMA(ddsShrink.11, alpha = 0.05)

# Volcano plots

volcanoTab.11 <- shrinkTab.11 %>%
  mutate(`-log10(pvalue)`= -log10(pvalue))

ggplot(volcanoTab.11, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = pvalue < 0.05), size = 1) +
  geom_text(data = ~top_n(.x, 1, wt=-FDR), aes(label=Symbol))

# Exercise 2

ddsShrink.33 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.33,
                          type = "ashr")
shrinkTab.33 <- as.data.frame(ddsShrink.33) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC = log2FoldChange, FDR = padj)

volcanoTab.33 <- shrinkTab.33 %>%
  mutate(`-log10(pvalue)`= -log10(pvalue))

ggplot(volcanoTab.33, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = pvalue < 0.05), size = 1) +
  geom_text(data = ~top_n(.x, 1, wt=-FDR), aes(label=Symbol))

# Venn Diagram

library(ggvenn)

vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>%
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange > 0) %>%
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange < 0) %>%
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange > 0) %>%
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange < 0)
vennDat

ggvenn(vennDat, set_name_size = 4)

# Heatmaps

library(ComplexHeatmap)
library(circlize)

sigGene <- shrinkTab.11 %>%
  top_n(300, wt=-FDR) %>%
  pull("GeneID")

plotDat <- vst(ddsObj.interaction)[sigGene,] %>%
  assay()

z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE)

ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")])

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        split = 3,
        rect_gp = gpar(col = "lightgrey", lwd = 0.3),
        top_annotation = ha1)

ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")], 
                        col = list(Status = c("Uninfected" = "hotpink", 
                                              "Infected" = "purple"), 
                                   TimePoint = c("d11" = "lightblue", 
                                                 "d33" = "darkblue")))
