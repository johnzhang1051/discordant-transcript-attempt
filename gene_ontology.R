# If not already installed
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

####### Get Data:

# Pull transcripts
transcripts <- read.csv("guide_effect/discordant_transcripts_non_overlaps.csv")

## extract to vector
symbols <- as.character(transcripts$Gene)

# Now run bitr with the vectors
entrez <- bitr(symbols, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)

## GO terms
gene_ontology_terms <- enrichGO(gene         = entrez$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9,
                    readable     = TRUE)

# View top results
head(gene_ontology_terms)

library(ggplot2)
dotplot(gene_ontology_terms, showCategory = 15) +
  ggtitle("GO Biological Processes Enriched")
barplot(gene_ontology_terms, showCategory = 15, title = "Top GO Terms")

## genes linked to specific GO
library(dplyr)
library(stringr)
# Search for the GO term (flexible search)
go_of_interest <- gene_ontology_terms %>%
  filter(str_detect(Description, "positive regulation of receptor clustering"))

# See genes linked
genes_linked <- strsplit(go_of_interest$geneID[1], "/")[[1]]

# View them
genes_linked


# Reactome for correlation_genes

## ignore signficiance
reactome_results <- enrichPathway(
  gene         = entrez$ENTREZID,
  organism     = "human",
  pvalueCutoff = 0.5, # include all results
  qvalueCutoff = 0.5, # include all results
  readable     = TRUE
)
# Convert to data frame
reactome_df <- as.data.frame(reactome_results)

dotplot(reactome_results, showCategory = 15) +
  ggtitle("Reactome Pathways Enriched")


## KEGG Pathway
# KEGG requires organism code — hsa = human
kegg_results <- enrichKEGG(gene = entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.5)

View(kegg_results@result)

dotplot(kegg_results, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")


# Filter for "MITF-M regulated melanocyte development"
mitf_pathway <- reactome_results@result %>%
  dplyr::filter(grepl("MITF", Description, ignore.case = TRUE))
