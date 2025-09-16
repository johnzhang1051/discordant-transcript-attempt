####################### Setting up environment

library(dplyr)

####################### Getting data
# Get transcript/isoform expression data - so we know which transcripts are expressed
transcript_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv")
)
colnames(transcript_expression_data)[1] <- "index"
# Filter out ACH-000931 and ACH-000008 because duplicates
duplicate_model_id <- c("ACH-000008", "ACH-000931")  # Add any others you found
transcript_expression_data <- transcript_expression_data %>%
  filter(!model_id %in% duplicate_model_id)

# Get protein 
protein_gene_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
)
colnames(protein_gene_expression_data)[1] <- "model_id"

####################### Filter models/cell-lines to melanoma and protein-encoding transcripts (do that later)
# Get cell line info - need this to identify melanoma cells
cell_info <- depmap_metadata()


# filter to all melanoma types that are not uveal, acral, mucosal, or uveal
specific_subtypes <- c("Melanoma", "Melanoma, amelanotic")
melanoma_cells <- cell_info %>%
  filter(subtype_disease %in% specific_subtypes)

# Subset for melanoma cell lines only
melanoma_transcript_expression <- transcript_expression_data %>%
  filter(model_id %in% melanoma_cells$depmap_id)
melanoma_gene_expression <- protein_gene_expression_data %>%
  filter(model_id %in% melanoma_cells$depmap_id)

# Subset gene expression to only samples in transcript data
melanoma_gene_expression <- melanoma_gene_expression %>%
  filter(model_id %in% melanoma_transcript_expression$model_id)

# FILTER TO PROTEIN ENCODING (TO DO)

write.csv(melanoma_transcript_expression, "cleaned_data/melanoma_transcript_expression.csv", row.names = FALSE)
write.csv(melanoma_gene_expression, "cleaned_data/melanoma_gene_expression.csv", row.names = FALSE)

####################### Identify transcript and gene associated with MITF-M expression

# Determine which col in transcript data is MITF-M "ENST00000394351"
MITF_transcript_column <- grep("ENST00000394351", colnames(transcript_expression_data), value = TRUE, ignore.case = TRUE)
print(MITF_transcript_column)


# Determine which col in gene data encodes for MITF-M
# using ensembl Biomart R package
# documentation: https://useast.ensembl.org/info/data/biomart/biomart_r_package.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)
# Function to automatically query biomaRt with attribute grouping

auto_biomart_query <- function(attributes, filters, values, mart) {
  all_attributes <- listAttributes(mart)
  attr_page_map <- all_attributes[all_attributes$name %in% attributes, c("name", "page")]
  pages <- split(attr_page_map$name, attr_page_map$page)
  
  results_list <- list()
  
  for (page_name in names(pages)) {
    page_attributes <- pages[[page_name]]
    
    if (!filters %in% page_attributes) {
      page_attributes <- c(filters, page_attributes)
    }
    
    tryCatch({
      page_result <- getBM(
        attributes = page_attributes,
        filters = filters,
        values = values,
        mart = mart
      )
      results_list[[page_name]] <- page_result
    }, error = function(e) {
      results_list[[page_name]] <- NULL
    })
  }
  
  return(results_list)
}

# Your attributes
attributes <- c(
  "ensembl_gene_id",
  "ensembl_gene_id_version", 
  "external_gene_name",
  "start_position",
  "end_position",
  "ensembl_transcript_id",
  "ensembl_transcript_id_version",
  "external_transcript_name"
)

# Remove version number from transcript ID - for example "ENST00000250838.6" -> "ENST00000250838"
transcript_ids <- sub("\\..*", "", sub("^[.]*", "", MITF_transcript_column))

ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Run the automated query
results <- auto_biomart_query(
  attributes = attributes,
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl_data
)

MITF_gene_column <- grep(results$feature_page$external_gene_name, colnames(protein_gene_expression_data), value = TRUE, ignore.case = TRUE)
print(MITF_gene_column)
