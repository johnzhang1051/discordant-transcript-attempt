####################### Setting up environment
library(depmap)
library(tidyverse)
library(dplyr)
library(GenomicRanges)

####################### Load transcript lists
discordant_transcripts <- readr::read_csv(
  file=file.path("correlation_results", "discordant_transcripts.csv")
)
colnames(discordant_transcripts)[1] <- "transcript_id"


correlated_transcripts <- readr::read_csv(
  file=file.path("correlation_results", "MITF_transcript_correlations_filtered.csv")
)


# Remove version numbers if present (e.g., "ENST00000394351.9" -> "ENST00000394351")
discordant_transcripts$transcript_id_clean <- sub("\\..*", "", discordant_transcripts$transcript_id)
correlated_transcripts$transcript_id_clean <- sub("\\..*", "", correlated_transcripts$transcript_id)

####### EDIT THIS
transcript_list <- discordant_transcripts
#transcript_list <- "correlated_transcripts"

transcript_list_name <- "discordant_transcripts"
#transcript_list_name <- "correlated_transcripts"



#######################  Query biomaRt for your transcripts
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")

library(biomaRt)

ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

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

# Your attributes - focusing on what's needed for guide analysis
attributes <- c(
  "ensembl_transcript_id",
  "external_gene_name",
  "ensembl_exon_id",
  "exon_chrom_start",
  "exon_chrom_end",
  "strand",
  "chromosome_name"
)

# Use your transcript list
transcript_ids <- unique(transcript_list$transcript_id_clean)

cat("Querying biomaRt for", length(transcript_ids), "transcripts...\n")

# Run the automated query
results <- auto_biomart_query(
  attributes = attributes,
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl_data
)

# export Biomart data
write.csv(results$structure, paste0("guide_effect/", transcript_list_name, "_exon_locations.csv"), row.names = FALSE)


#######################  Get exon-level dataframe
exon_locations <- results$structure

exon_locations <- readr::read_csv(
  file=paste0("guide_effect/", transcript_list_name, "_exon_locations.csv")
)

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
  dplyr::select(-strand)  # Remove conflicting column names

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
  dplyr::select(-guide_strand)  # Remove original strand column

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
write.csv(overlaps, paste0("guide_effect/", transcript_list_name, "_exon_guide_overlaps.csv"), row.names = FALSE)


###################### Look at non-overlaps

transcripts_with_guides <- unique(overlaps$exon_data.ensembl_transcript_id)
non_overlaps <- transcript_list %>%
  filter(!transcript_id_clean %in% transcripts_with_guides)

write.csv(non_overlaps, paste0("guide_effect/", transcript_list_name, "_non_overlaps.csv"), row.names = FALSE)


#######################  Correlate Guide Effect to Transcripts

# Pull in guide effect data from Depmap
avana_raw_readcounts <- readr::read_csv(
  file=file.path("depmap-data", file="AvanaRawReadcounts.csv")
)
colnames(avana_raw_readcounts)[1] <- "sg_rna"

screen_sequence_map <- readr::read_csv(
  file=file.path("depmap-data", file="ScreenSequenceMap.csv")
)

avana_effects <- avana_raw_readcounts %>%
  pivot_longer(cols = -sg_rna, names_to = "SequenceID", values_to = "read_count") %>%
  # Calculate guide depletion (lower counts = more lethal)
  group_by(SequenceID) %>%
  mutate(
    # Log2 transform and center by median (common pre-processing)
    log_count = log2(read_count + 1),
    guide_effect = log_count - median(log_count, na.rm = TRUE)
  ) %>%
  ungroup()

# aggregate guide effects to only hit melanoma types

# Get cell line info - need this to identify melanoma cells
cell_info <- depmap_metadata()

# filter to all melanoma types that are not uveal, acral, mucosal, or uveal
specific_subtypes <- c("Melanoma", "Melanoma, amelanotic")
melanoma_cells <- cell_info %>%
  filter(subtype_disease %in% specific_subtypes)

melanoma_cells_sequences <- melanoma_cells %>% 
  left_join(screen_sequence_map, by = c("depmap_id" = "ModelID"))

melanoma_cells_sequences <- unique(melanoma_cells_sequences %>% dplyr::select(SequenceID, subtype_disease))

melanoma_effects <- avana_effects %>% 
  inner_join(melanoma_cells_sequences, by = c("SequenceID" = "SequenceID"))

melanoma_effects_aggregated <- melanoma_effects %>%
  group_by(sg_rna) %>%
  summarise(
    mean_read_count = mean(read_count, na.rm = TRUE),
    mean_log_count = mean(log_count, na.rm = TRUE),
    mean_guide_effect = mean(guide_effect, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

transcript_guide_effects <- overlaps %>% 
  left_join(melanoma_effects_aggregated, by = c("guide_data.sgRNA" = "sg_rna"))

# average guide effects by transcripts_id
transcript_aggregated_effects <- transcript_guide_effects %>%
  group_by(exon_data.ensembl_transcript_id) %>%
  summarise(
    gene_name = dplyr::first(exon_data.external_gene_name),
    n_guides = n(),
    mean_guide_effect = mean(mean_guide_effect, na.rm = TRUE),
    median_guide_effect = median(mean_guide_effect, na.rm = TRUE),
    sd_guide_effect = sd(mean_guide_effect, na.rm = TRUE),
    min_guide_effect = min(mean_guide_effect, na.rm = TRUE),
    max_guide_effect = max(mean_guide_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_guide_effect))

# export final results
write.csv(transcript_aggregated_effects,paste0("guide_effect/", transcript_list_name, "_guide_effects.csv"), row.names = FALSE)


library(ggplot2)

# Create boxplot of guide effects
guide_effects_plot <- ggplot(transcript_aggregated_effects, aes(x = "All Transcripts", y = mean_guide_effect)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(
    title = paste("Guide Effects for", str_to_title(transcript_list_name), "Transcripts"),
    x = "",
    y = "Mean Guide Effect",
    subtitle = paste("n =", nrow(transcript_aggregated_effects), "transcripts")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)

# Display the plot
print(guide_effects_plot)
