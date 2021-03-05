library(biomaRt)
library(tidyverse)

results.interaction.11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")

## set up connection to ensembl database
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'mmusculus_gene_ensembl',
                      version = 102)


# get annot
filterType <- "ensembl_gene_id"
# get the Ensembl IDs from our results table
filterValues <- rownames(results.interaction.11)

# run the query
attributeNames <- c("ensembl_gene_id",
                    "entrezgene_id",
                    "external_gene_name",
                    "description",
                    "gene_biotype",
                    "chromosome_name",
                    "start_position",
                    "end_position",
                    "strand",
                    "entrezgene_accession")

# Get annotations
annot <- getBM(attributes=attributeNames,
               filters = filterType,
               values = filterValues,
               mart = ensembl)

saveRDS(annot, file="Full_annotation_with_duplicates.rds")

#annot <- readRDS(file="Full_annotation_with_duplicates.rds")

# There are ensembl id's with multiple Entrez ID's
# Deduplicate the entrez IDS by matching the entrez symbol and the 
# "external_gene_name"

dups <- annot %>%  
    add_count(ensembl_gene_id) %>%  
    filter(n>1)

fixedDuplicates <- dups %>% 
    select(-n) %>% 
    filter(entrezgene_accession==external_gene_name)

# check that this has no duplicates
fixedDuplicates %>%  
    add_count(ensembl_gene_id) %>%  
    filter(n>1)
## NONE

annot2 <- annot %>%  
    add_count(ensembl_gene_id) %>%  
    filter(n==1) %>% 
    select(-n) %>% 
    bind_rows(fixedDuplicates)

nrow(annot2)
length(unique(annot$ensembl_gene_id))

# for the four remaining just arbitrary decsion
fixedDuplicates <- dups %>% 
    filter(!ensembl_gene_id%in%annot2$ensembl_gene_id) %>% 
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    select(-n)

annotUn <- bind_rows(annot2, fixedDuplicates)
nrow(annotUn)
length(unique(annot$ensembl_gene_id))
all(filterValues%in%annotUn$ensembl_gene_id)

# The problem: we now have 31 duplicated Entrez IDs
annotUn %>% 
    filter(!is.na(entrezgene_id)) %>% 
    add_count(entrezgene_id) %>% 
    filter(n>1) %>% 
    dplyr::count(entrezgene_id)

# fsgea throws a nasty warning about multiple genes if we have this issue
as.data.frame(results.interaction.11) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(annotUn) %>% 
    filter(!is.na(entrezgene_id)) %>% 
    add_count(entrezgene_id) %>% 
    filter(n>1) %>% 
    select(ensembl_gene_id, entrezgene_id, padj, chromosome_name) %>% 
    arrange(padj)

# we need a pragmatic solution for the course
# these genes are mostly non-significant, we'll arbitrarily set the
# second entry to NA. Some of the duplicates are on patch scaffolds, we'll 
# arrange by chromosome, so that these get set to NA
dupEntrez <- annotUn %>% 
    add_count(entrezgene_id) %>% 
    filter(n>1) %>% 
    select(-n) %>% 
    arrange(entrezgene_id, chromosome_name)
dupEntrez$entrezgene_id[duplicated(dupEntrez$entrezgene_id)] <- NA    


annotFin <- annotUn %>% 
    add_count(entrezgene_id) %>% 
    filter(n==1) %>% 
    select(-n) %>% 
    bind_rows(dupEntrez)

dim(annotFin)
annotFin %>% 
    filter(!is.na(entrezgene_id)) %>% 
    add_count(entrezgene_id) %>% 
    filter(n>1) 
all(filterValues%in%annotFin$ensembl_gene_id)

### Final table

ensemblAnnot <- rownames(results.interaction.11) %>%  
    enframe(name = NULL, value = "ensembl_gene_id")  %>%  
    left_join(annotFin) %>%
    dplyr::select(GeneID="ensembl_gene_id", Entrez="entrezgene_id",
                  Symbol="external_gene_name", Description="description",
                  Biotype="gene_biotype", Chr="chromosome_name",
                  Start="start_position", End="end_position",
                  Strand="strand")

saveRDS(ensemblAnnot, file="RObjects/Ensembl_annotations.rds")
