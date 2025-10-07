# If not already installed
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(dplyr)
library(stringr)

####### Get Data:

# Pull transcripts
transcripts <- read.csv("guide_effect/mitf_high/discordant_MITF_HIGH_guide_effect_results.csv")

# filter to non-screened transcripts
transcripts <- transcripts %>% filter(n_exons_targeted == 0) 

## extract to vector
symbols <- as.character(transcripts$gene_name)

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

go_plot <- dotplot(gene_ontology_terms, showCategory = 15) +
  ggtitle("GO Biological Processes Enriched")

ggsave("gene_ontology/go_process.png", go_plot, width = 10, height = 8, dpi = 300)

## genes linked to specific GO
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

reactome_plot <- dotplot(reactome_results, showCategory = 15) +
  ggtitle("Reactome Pathways Enriched")

ggsave("gene_ontology/reactome_pathways.png", reactome_plot, width = 10, height = 8, dpi = 300)


## KEGG Pathway
# KEGG requires organism code — hsa = human
kegg_results <- enrichKEGG(gene = entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.5)

View(kegg_results@result)

kegg_plot <- dotplot(kegg_results, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")

ggsave("gene_ontology/kegg_pathways.png", kegg_plot, width = 10, height = 8, dpi = 300)

# Filter for "MITF-M regulated melanocyte development"
mitf_pathway <- reactome_results@result %>%
  dplyr::filter(grepl("MITF", Description, ignore.case = TRUE))
