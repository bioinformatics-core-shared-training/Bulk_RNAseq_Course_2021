library(DESeq2)
library(tidyverse)

# Load the data

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv")

colnames(txi$counts)==sampleinfo$SampleName

# Create DESeq2 object

## create model

# y = 2x
# y ~ 2x

simple.model <- as.formula( ~ Status)

model.matrix(simple.model, data = sampleinfo)

sampleinfo <- sampleinfo %>% 
  mutate(Status = fct_relevel(Status, "Uninfected"))

model.matrix(simple.model, data = sampleinfo)

## Create DESeqDataSet

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)

## Filter genes

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

nrow(ddsObj.raw)
nrow(ddsObj.filt)

# DESeq2 Workflow


# estimate the size factors

ddsObj <- estimateSizeFactors(ddsObj.filt) 

normalizationFactors(ddsObj.filt)

head(normalizationFactors(ddsObj))

## MA plot

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)

limma::plotMA(logcounts, array = 5, ylim = c(-5, 5))
abline(h = 0, col="red")


logcounts <- log2(counts(ddsObj, normalized = TRUE) + 1)

limma::plotMA(logcounts, array = 5, ylim = c(-5, 5))
abline(h = 0, col="red")


# estimate dispersion factors

ddsObj <- estimateDispersions(ddsObj)


plotDispEsts(ddsObj)

# Apply the Wald test

ddsObj <- nbinomWaldTest(ddsObj)

# The DESeq command

ddsObj <- DESeq(ddsObj.filt)

# Generate a results table 

results.simple <- results(ddsObj, alpha = 0.05)
results.simple

# How many genes are differentially expressed

sum(results.simple$padj < 0.05)
View(as.data.frame(results.simple))


sum(is.na(results.simple$padj))


sum(results.simple$padj < 0.05, na.rm = TRUE)


### upregulated

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE)

### downregulated

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE)


# additive model

additive.model <- as.formula(~ TimePoint + Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# run the DESeq2 workflow

ddsObj <- DESeq(ddsObj.filt)

# extract the results

results.additive <- results(ddsObj, alpha = 0.05)

results.additive

sum(results.additive$padj < 0.05, na.rm = TRUE)


# The default contrast

model.matrix(additive.model, data=sampleinfo)

resultsNames(ddsObj)

results.InfectedvUninfected <- results.additive
rm(results.additive)

top100genesIvU <- as.data.frame(results.InfectedvUninfected) %>% 
  rownames_to_column("GeneID") %>% 
  top_n(100, wt = -padj)

# Exercise 2

resultsNames(ddsObj)

results.d33vd11 <- results(ddsObj, 
                           name = "TimePoint_d33_vs_d11",
                           alpha = 0.05)

results.d33vd11

sum(results.d33vd11$padj < 0.05, na.rm=TRUE)



# PCA

vstcounts <- vst(ddsObj, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))


# Comparing two design models

# additive v simple

ddsObj.LRT <- DESeq(ddsObj, test = "LRT", reduced = simple.model)
results.Additive_v_Simple <- results(ddsObj.LRT)

results.Additive_v_Simple

sum(results.Additive_v_Simple$padj < 0.05, na.rm=TRUE)
# 66

# interaction v additive

interaction.model <- as.formula(~ TimePoint + Status + TimePoint:Status)
interaction.model <- as.formula(~ TimePoint * Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.interaction <- DESeq(ddsObj.filt)

ddsObj.LRT <- DESeq(ddsObj.interaction, test="LRT", reduced=additive.model)

results.InteractionvAdditive <- results(ddsObj.LRT)
table(results.InteractionvAdditive$padj < 0.05)


# Extracting specific contrasts from the interactions model


resultsNames(ddsObj.interaction)

results.interaction.d11 <- results(ddsObj.interaction,
                                   name = "Status_Infected_vs_Uninfected",
                                   alpha = 0.05)

results.interaction.d11 


## day 33 contrast

results.interaction.d33 <- results(ddsObj.interaction,
                                   contrast =list(c("Status_Infected_vs_Uninfected",
                                                    "TimePointd33.StatusInfected")),
                                   alpha = 0.05)
results.interaction.d33


# how many genes DE at d11

sum(results.interaction.d11$padj < 0.05, na.rm = TRUE)
# 1072

# how many genes DE at d33

sum(results.interaction.d33$padj < 0.05, na.rm = TRUE)
# 2782

# Exercise 4

# uninfected mice

resultsNames(ddsObj.interaction)

results.d33vd11_uninfected <- results(ddsObj.interaction,
                                      name = "TimePoint_d33_vs_d11",
                                      alpha = 0.05)
table(results.d33vd11_uninfected$padj < 0.05)

# infected mice


resultsNames(ddsObj.interaction)

results.d33vd11_infected <- results(ddsObj.interaction,
                                      contrast = list(c("TimePoint_d33_vs_d11",
                                                        "TimePointd33.StatusInfected")),
                                      alpha = 0.05)
table(results.d33vd11_infected$padj < 0.05)







