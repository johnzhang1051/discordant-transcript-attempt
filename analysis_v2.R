####################### Installing depmap

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("depmap")

####################### Setting up environment

library(depmap)
library(tidyverse)
library(dplyr)
library(corrr)
library(ggplot2)
library(tidyr)
setwd("/Users/johnz/Documents/GitFiles/discordant-transcript-attempt")

####################### Getting data

# 1. Get cell line info - need this to identify melanoma cells (actually don't need this)
cell_info <- depmap_metadata()
melanoma_cells <- filter(cell_info, grepl("Melanoma", subtype_disease, ignore.case = TRUE))

# 2. Get transcript/isoform expression data - so we know which transcripts are expressed
transcript_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv")
)
colnames(transcript_expression_data)[1] <- "index"

# 3. Get protein 
protein_gene_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
)
colnames(protein_gene_expression_data)[1] <- "model_id"
protein_gene_expression_data <- protein_gene_expression_data %>%
  rename("MITF" = "MITF (4286)")

####################### Filter data

# Get MITF expression data:
MITF_columns <- grep("MITF", colnames(protein_gene_expression_data), value = TRUE, ignore.case = TRUE)
MITF_expression_data <- protein_gene_expression_data %>%
  select(model_id, all_of(MITF_columns))

####################### Combine data
combined_data <- MITF_expression_data %>%
  inner_join(transcript_expression_data, by = "model_id")

####################### Correlation analysis between MITF and transcripts
transcript_names <- colnames(combined_data)[!colnames(combined_data) %in% c("model_id", "MITF", "index", "profile_id", "is_default_entry")]

total_transcripts <- length(transcript_names)

MITF_transcript_correlations <- map_dfr(transcript_names, function(transcript) {
  # Get current position
  current_pos <- which(transcript_names == transcript)
  
  # Print progress every 100 transcripts
  if(current_pos %% 100 == 0) {
    cat("Correlations finished/total:", current_pos, "/", total_transcripts, 
        "(", round(current_pos/total_transcripts*100, 1), "%)\n")
  }
  
  tibble(
    term = transcript,
    pearson_corr = cor(combined_data$MITF, combined_data[[transcript]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(combined_data$MITF, combined_data[[transcript]], 
                        method = "spearman", use = "complete.obs")
  )
})

# export correlations because they take forever to run
write.csv(MITF_transcript_correlations,"MITF_transcript_correlations.csv", row.names = FALSE)

####################### Correlation analysis for gene-level expression
transcript_names <- colnames(protein_gene_expression_data)[!colnames(protein_gene_expression_data) %in% c("model_id")]

total_transcripts <- length(transcript_names)

MITF_gene_correlations <- map_dfr(transcript_names, function(transcript) {
  # Get current position
  current_pos <- which(transcript_names == transcript)
  
  # Print progress every 100 transcripts
  if(current_pos %% 100 == 0) {
    cat("Correlations finished/total:", current_pos, "/", total_transcripts, 
        "(", round(current_pos/total_transcripts*100, 1), "%)\n")
  }
  
  # tibble creates a small DF with the results (in this case the columns are: term, pearson_corr, spearman_corr)
  tibble(
    term = transcript,
    # Run correlations between MITF<->a transcript
    # complete.obs: removes pairs where any data is missing
    pearson_corr = cor(protein_gene_expression_data$MITF, protein_gene_expression_data[[transcript]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(protein_gene_expression_data$MITF, protein_gene_expression_data[[transcript]], 
                        method = "spearman", use = "complete.obs")
  )
})

# export correlations because they take forever to run
write.csv(MITF_gene_correlations,"MITF_gene_correlations.csv", row.names = FALSE)


#######################  Pull in correlations data (if working on it)
MITF_transcript_correlations <- readr::read_csv(
  file="MITF_transcript_correlations.csv"
)

colnames(MITF_transcript_correlations)[1] <- "transcript_id"

MITF_gene_correlations <- readr::read_csv(
  file="MITF_gene_correlations.csv"
)
colnames(MITF_gene_correlations)[1] <- "gene_id"

filtered_transcripts <- MITF_transcript_correlations %>%
  filter(pearson_corr > 0 | spearman_corr > 0)

#######################  Map transcripts to parent genes (not sure how to do this yet)

# using ensembl Biomart R package
# documentation: https://useast.ensembl.org/info/data/biomart/biomart_r_package.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

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
  "external_transcript_name",
  "transcript_start",
  "transcript_end", 
  "ensembl_exon_id",
  "hgnc_symbol",
  "exon_chrom_start",
  "exon_chrom_end",
  "strand",
  "chromosome_name"
)

# Remove version number from transcript ID - for example "ENST00000250838.6" -> "ENST00000250838"
transcript_ids <- sub("\\..*", "", sub("^[.]*", "", filtered_transcripts$transcript_id))

# Run the automated query
results <- auto_biomart_query(
  attributes = attributes,
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl_data
)

# export Biomart data because they take forever to run
write.csv(results$feature_page,"biomart_feature_page.csv", row.names = FALSE)
write.csv(results$homologs,"biomart_homologs.csv", row.names = FALSE)
write.csv(results$sequences,"biomart_sequences.csv", row.names = FALSE)
write.csv(results$snp,"biomart_snp.csv", row.names = FALSE)
write.csv(results$snp_somatic,"biomart_snp_somatic.csv", row.names = FALSE)
write.csv(results$structure,"biomart_structure.csv", row.names = FALSE)

#######################  Get exon-level dataframe
exon_locations <- readr::read_csv(
  file="biomart_structure.csv"
)
#exon_locations <- results$structure

#######################  Map Depmap Guides to Transcripts

# Use Avana data check if sgRNAs overlap with Exon coordinates 
avana_guide_map <- readr::read_csv(
  file=file.path("depmap-data", file="AvanaGuideMap.csv")
)

# filter to coordinates without "NT", and only have "chr..."
avana_guide_map <- avana_guide_map %>%
  filter(grepl("^chr[0-9XY]", GenomeAlignment))

# extract sgRNA targeted coordinate and strand
avana_guide_map <- avana_guide_map %>%
  separate(GenomeAlignment, 
           into = c("chromosome", "guide_coordinate", "guide_strand"), 
           sep = "_",
           remove = FALSE) %>%
  mutate(
    chromosome = gsub("chr", "", chromosome),
    guide_coordinate = as.integer(guide_coordinate),
    guide_strand = case_when(
      guide_strand == "+" ~ 1,
      guide_strand == "-" ~ -1,
      TRUE ~ NA_real_
    )
  )


# Use GenomicRanges library to find overlaps between exon and guide coordinates for same chromosomes
library(GenomicRanges)

# Fix exon data by removing conflicting columns and converting strand
exon_locations <- exon_locations %>%
  mutate(strand_char = case_when(
    strand == 1 ~ "+",
    strand == -1 ~ "-",
    TRUE ~ "*"
  )) %>%
  select(-strand)  # Remove conflicting column names

exon_gr <- makeGRangesFromDataFrame(
  exon_locations,
  seqnames.field = "chromosome_name",
  start.field = "exon_chrom_start", 
  end.field = "exon_chrom_end",
  strand.field = "strand_char",
  keep.extra.columns = TRUE
)

# Fix guide data similarly
avana_guide_map <- avana_guide_map %>%
  mutate(strand_char = case_when(
    guide_strand == 1 ~ "+",
    guide_strand == -1 ~ "-", 
    TRUE ~ "*"
  )) %>%
  select(-guide_strand)  # Remove original strand column, because GenomicRange dhas it's own columns it needs

guide_gr <- makeGRangesFromDataFrame(
  avana_guide_map,
  seqnames.field = "chromosome",
  start.field = "guide_coordinate",
  end.field = "guide_coordinate",
  strand.field = "strand_char",
  keep.extra.columns = TRUE
)

# Find overlaps
overlap_hits <- findOverlaps(guide_gr, exon_gr)

# Extract overlapping data
overlaps <- data.frame(
  guide_data = avana_guide_map[queryHits(overlap_hits), ],
  exon_data = exon_locations[subjectHits(overlap_hits), ]
)


# export overlaps data because it takes a while to run
write.csv(overlaps,"exon_guide_overlaps.csv", row.names = FALSE)

#######################  Get overlap data (START HERE JOHN)
overlaps <- readr::read_csv(
  file="exon_guide_overlaps.csv"
)

#######################  Correlate Guide Effect to Transcripts

# Pull in guide effect data from Depmap


