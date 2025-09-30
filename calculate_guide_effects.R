####################### Setting up environment
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

####################### CONFIGURABLE INPUT - CHOOSE YOUR TRANSCRIPT LIST
# Modify these two lines to switch between different transcript lists:
transcript_list_file <- "correlation_results/discordant_transcripts.csv"  # Changed to discordant
transcript_list_name <- "discordant_transcripts"  # Changed to discordant

# Load the chosen transcript list
transcript_list <- readr::read_csv(file = transcript_list_file)
colnames(transcript_list)[1] <- "transcript_id"

# Remove version numbers if present
transcript_list$transcript_id_clean <- sub("\\..*", "", transcript_list$transcript_id)

cat("Analyzing", nrow(transcript_list), transcript_list_name, "transcripts\n")

#######################  Query biomaRt for your transcripts
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
overlap_hits <- findOverlaps(guide_gr, exon_gr, ignore.strand = TRUE)  # Added ignore.strand

# Extract overlapping data
overlaps <- data.frame(
  guide_data = avana_guide_map[queryHits(overlap_hits), ],
  exon_data = exon_locations[subjectHits(overlap_hits), ]
)

# export overlaps data because it takes a while to run
write.csv(overlaps, paste0("guide_effect/", transcript_list_name, "_exon_guide_overlaps.csv"), row.names = FALSE)

###################### TRANSCRIPT-LEVEL GUIDE COUNTS

# Count how many guides target each transcript
transcript_guide_counts <- overlaps %>%
  group_by(exon_data.ensembl_transcript_id) %>%
  summarise(
    gene_name = dplyr::first(exon_data.external_gene_name),
    n_guides_targeting = n_distinct(guide_data.sgRNA),
    n_exons_targeted = n_distinct(exon_data.ensembl_exon_id),
    guide_list = paste(unique(guide_data.sgRNA), collapse = ";"),
    .groups = "drop"
  ) %>%
  rename(transcript_id = exon_data.ensembl_transcript_id)

# Get transcripts with guides vs without guides
transcripts_with_guides <- unique(overlaps$exon_data.ensembl_transcript_id)
transcripts_without_guides <- transcript_list %>%
  filter(!transcript_id_clean %in% transcripts_with_guides) %>%
  select(transcript_id = transcript_id_clean) %>%
  mutate(
    gene_name = NA,
    n_guides_targeting = 0,
    n_exons_targeted = 0,
    guide_list = ""
  )

# Combine to get complete transcript guide summary
transcript_guide_summary <- bind_rows(
  transcript_guide_counts,
  transcripts_without_guides
)

# Export transcript-level guide counts
write.csv(transcript_guide_summary, 
          paste0("guide_effect/", transcript_list_name, "_transcript_guide_counts.csv"), 
          row.names = FALSE)

###################### GUIDE-LEVEL ANALYSIS

# Create guide-level dataframe showing which transcripts each guide targets
guide_level_analysis <- overlaps %>%
  group_by(guide_data.sgRNA) %>%
  summarise(
    n_transcripts_targeted = n_distinct(exon_data.ensembl_transcript_id),
    n_genes_targeted = n_distinct(exon_data.external_gene_name),
    n_exons_hit = n(),
    transcript_list = paste(unique(exon_data.ensembl_transcript_id), collapse = ";"),
    gene_list = paste(unique(exon_data.external_gene_name), collapse = ";"),
    chromosomes = paste(unique(exon_data.chromosome_name), collapse = ";"),
    .groups = "drop"
  ) %>%
  rename(sgRNA = guide_data.sgRNA)

# Export guide-level analysis
write.csv(guide_level_analysis, 
          paste0("guide_effect/", transcript_list_name, "_guide_level_analysis.csv"), 
          row.names = FALSE)

###################### Look at non-overlaps
non_overlaps <- transcript_list %>%
  filter(!transcript_id_clean %in% transcripts_with_guides)

write.csv(non_overlaps, paste0("guide_effect/", transcript_list_name, "_non_overlaps.csv"), row.names = FALSE)

#######################  Correlate Guide Effect to Transcripts - SPLIT BY MELANOMA VS NON-MELANOMA

# Pull in guide effect data from Depmap
avana_logfold_change <- readr::read_csv(
  file=file.path("depmap-data", file="AvanaLogfoldChange.csv")
)
colnames(avana_logfold_change)[1] <- "sgRNA"

screen_sequence_map <- readr::read_csv(
  file=file.path("depmap-data", file="ScreenSequenceMap.csv")
)

# Get cell line info
cell_info <- read.csv("depmap-data/Model.csv")

# Define melanoma types
melanoma_subtypes <- c("Melanoma", "Cutaneous Melanoma")

# Filter to only guides that target our transcripts FIRST
guides_of_interest <- unique(overlaps$guide_data.sgRNA)


# Only pivot the guides we care about
filtered_lfc <- avana_logfold_change %>%
  filter(sgRNA %in% guides_of_interest)

# Now pivot - much smaller dataset
guide_effects <- filtered_lfc %>%
  pivot_longer(cols = -sgRNA, 
               names_to = "SequenceID", 
               values_to = "log_fold_change") %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID") %>%
  filter(!is.na(OncotreeSubtype)) %>%  # Remove rows without cancer type
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

# Join with transcript overlaps
transcript_guide_effects <- overlaps %>%
  left_join(
    guide_effects,
    by = c("guide_data.sgRNA" = "sgRNA"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(log_fold_change)) %>%
  # Remove duplicate guide-transcript-cell line combinations
  distinct(exon_data.ensembl_transcript_id, guide_data.sgRNA, SequenceID, .keep_all = TRUE)

# Aggregate by transcript and cancer type
transcript_effects_by_cancer <- transcript_guide_effects %>%
  group_by(exon_data.ensembl_transcript_id, is_melanoma) %>%
  summarise(
    gene_name = first(exon_data.external_gene_name),
    n_guides = n_distinct(guide_data.sgRNA),
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    median_lfc = median(log_fold_change, na.rm = TRUE),
    sd_lfc = sd(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

# Separate melanoma vs non-melanoma
melanoma_effects <- transcript_effects_by_cancer %>%
  filter(is_melanoma) %>%
  select(transcript_id = exon_data.ensembl_transcript_id,
         gene_name,
         n_guides,
         guide_effect_melanoma = mean_lfc,
         melanoma_median = median_lfc,
         melanoma_n_cell_lines = n_cell_lines)

non_melanoma_effects <- transcript_effects_by_cancer %>%
  filter(!is_melanoma) %>%
  select(transcript_id = exon_data.ensembl_transcript_id,
         guide_effect_non_melanoma = mean_lfc,
         non_melanoma_median = median_lfc,
         non_melanoma_n_cell_lines = n_cell_lines)

# Combine
transcript_comparison <- melanoma_effects %>%
  left_join(non_melanoma_effects, by = "transcript_id") %>%
  mutate(
    difference = guide_effect_melanoma - guide_effect_non_melanoma
  ) %>%
  arrange(difference)

###### Try to score the differences to identify interesting transcripts

transcript_comparison <- transcript_comparison %>%
  mutate(
    # Z-score style - how many SDs is melanoma from non-melanoma
    # Negative = melanoma-specific
    melanoma_specificity_zscore = (guide_effect_melanoma - guide_effect_non_melanoma) / 
      sqrt(melanoma_median^2 + non_melanoma_median^2 + 0.01),
    
    # Selectivity index (drug-like metric)
    # Positive = melanoma-specific
    selectivity_index = (guide_effect_non_melanoma - guide_effect_melanoma) / 
      (abs(guide_effect_non_melanoma) + abs(guide_effect_melanoma) + 0.01),
  )

# Z-Score < 1st quartile
# Selectivity > 3rd quartile
# Guide Effect < Median
zscore_q1 <- quantile(transcript_comparison$melanoma_specificity_zscore, 0.25, na.rm = TRUE)
selectivity_q3 <- quantile(transcript_comparison$selectivity_index, 0.75, na.rm = TRUE)
guide_effect_median <- median(transcript_comparison$guide_effect_melanoma, na.rm = TRUE)

# Identify interesting transcripts
interesting_transcripts <- transcript_comparison %>%
  filter(
    melanoma_specificity_zscore < zscore_q1  # Most negative z-scores
    & selectivity_index > selectivity_q3          # Highest selectivity
    & guide_effect_melanoma < guide_effect_median          # Highest selectivity
  ) %>%
  arrange(desc(selectivity_index)) %>% 
  select(transcript_id, gene_name, n_guides, 
                                          guide_effect_melanoma, guide_effect_non_melanoma,
                                          selectivity_index, melanoma_specificity_zscore)

View(interesting_transcripts)

# Add a flag to transcript_comparison for plotting
transcript_comparison <- transcript_comparison %>%
  mutate(is_interesting = transcript_id %in% interesting_transcripts$transcript_id)

# Export final results
write.csv(transcript_comparison, 
          paste0("guide_effect/", transcript_list_name, "_melanoma_vs_nonmelanoma.csv"), 
          row.names = FALSE)

###################### Graphs

# Scatter plot with interesting transcripts highlighted
ggplot(transcript_comparison,
       aes(x = guide_effect_non_melanoma, 
           y = guide_effect_melanoma,
           color = is_interesting,
           size = n_guides)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
  ggrepel::geom_text_repel(
    data = filter(transcript_comparison, is_interesting),
    aes(label = gene_name),
    size = 3,
    max.overlaps = 15,
    color = "darkred"
  ) +
  labs(
    title = paste(str_to_title(transcript_list_name), ": Melanoma vs Non-Melanoma"),
    subtitle = "Red points = melanoma-specific targets (below diagonal)",
    x = "Guide Effect in Non-Melanoma",
    y = "Guide Effect in Melanoma",
    color = "Interesting",
    size = "Number of Guides"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70"),
                     labels = c("TRUE" = "Melanoma-specific", "FALSE" = "Other"))

# 3. Difference plot with interesting transcripts highlighted
ggplot(transcript_comparison, 
       aes(x = reorder(transcript_id, difference), 
           y = difference,
           fill = is_interesting)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Melanoma Specificity",
    subtitle = "Red = interesting targets; Negative = guides more lethal in melanoma",
    x = "Transcript ID",
    y = "Difference (Melanoma - Non-Melanoma)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "Interesting", "FALSE" = "Other")) +
  coord_flip()
