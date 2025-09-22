####################### Setting up environment
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)
library(GenomicRanges)

# Fix select() conflicts
conflict_prefer("select", "dplyr")
conflicts_prefer(GenomicRanges::setdiff)
conflict_prefer("filter", "dplyr")

####################### Load Required Data

# Load the complete discordant transcript list (before filtering for overlaps)
all_discordant_transcripts <- readr::read_csv(
  file = file.path("correlation_results", "discordant_transcripts.csv")
)
colnames(all_discordant_transcripts)[1] <- "transcript_id"

# Load the non-overlapping transcripts (those NOT targeted by CRISPR)
non_overlap_transcripts <- readr::read_csv(
  file = file.path("guide_effect", "discordant_transcripts_non_overlaps.csv")
)
colnames(non_overlap_transcripts)[1] <- "transcript_id"

# Load overlaps data (transcript-guide mappings - these are the ones WITH CRISPR coverage)
overlaps <- readr::read_csv(
  file = "guide_effect/discordant_transcripts_exon_guide_overlaps.csv"
)

# Get the transcripts that DO have CRISPR coverage
transcripts_with_crispr <- unique(overlaps$exon_data.ensembl_transcript_id)

####################### Load ChIP-seq data for unique promoter analysis
Kenny <- read.csv("cleaned_data/Kenny.csv")
Laurette <- read.csv("cleaned_data/Laurette.csv") 
Louph <- read.csv("cleaned_data/Louphrasitthiphol.csv")

# Standardize transcript IDs in ChIP-seq data
standardize_transcripts <- function(df) {
  df %>%
    dplyr::rename(transcript_id = transcriptId) %>%
    dplyr::mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))
}

Kenny <- standardize_transcripts(Kenny)
Laurette <- standardize_transcripts(Laurette)
Louph <- standardize_transcripts(Louph)

####################### UNIQUE PROMOTER ANALYSIS

# Combine coordinates from all three datasets
all_transcript_coords <- bind_rows(
  Kenny %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Kenny"),
  Laurette %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Laurette"), 
  Louph %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Louph")
) %>%
  distinct(transcript_id, .keep_all = TRUE)  # Keep first occurrence if transcript appears in multiple datasets

# Extract transcript coordinates from Kenny dataset for promoter analysis
transcript_coords <- all_transcript_coords %>%
  dplyr::select(transcript_id, seqnames, start, end, strand, geneId) %>%
  distinct() %>%
  mutate(
    # Define promoter as 1kb upstream of TSS
    promoter_start = ifelse(strand == 1, start - 1000, end),
    promoter_end = ifelse(strand == 1, start, end + 1000),
    promoter_start = pmax(1, promoter_start)
  )

transcript_coords_clean <- transcript_coords %>%
  dplyr::select(-start, -end)

# Create GenomicRanges object for promoter regions
promoter_gr <- makeGRangesFromDataFrame(
  transcript_coords_clean,
  seqnames.field = "seqnames",
  start.field = "promoter_start", 
  end.field = "promoter_end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

# Find overlaps between promoter regions
overlaps_promoter <- findOverlaps(promoter_gr, promoter_gr)
overlap_summary <- as.data.frame(overlaps_promoter) %>%
  filter(queryHits != subjectHits)  # Exclude self-overlaps

# Get all transcripts that have ANY overlap with other transcripts
transcripts_with_overlaps <- unique(c(overlap_summary$queryHits, overlap_summary$subjectHits))

# Identify transcripts with unique (non-overlapping) promoters
unique_promoter_indices <- setdiff(1:length(promoter_gr), transcripts_with_overlaps)
unique_promoter_transcripts <- transcript_coords$transcript_id[unique_promoter_indices]


####################### ANALYZE NON-OVERLAPPING DISCORDANT TRANSCRIPTS

# Analyze non-overlapping transcripts for unique promoters
non_overlap_with_promoter_info <- non_overlap_transcripts %>%
  mutate(
    has_unique_promoter = transcript_id %in% unique_promoter_transcripts,
    has_coordinate_data = transcript_id %in% transcript_coords$transcript_id,
    has_crispr_coverage = transcript_id %in% transcripts_with_crispr
  )

# Summary statistics for non-overlapping transcripts
non_overlap_summary <- non_overlap_with_promoter_info %>%
  summarise(
    total_non_overlap = n(),
    with_coordinate_data = sum(has_coordinate_data, na.rm = TRUE),
    with_unique_promoter = sum(has_unique_promoter, na.rm = TRUE),
    with_overlapping_promoter = sum(!has_unique_promoter & has_coordinate_data, na.rm = TRUE),
    without_coordinate_data = sum(!has_coordinate_data, na.rm = TRUE),
    proportion_unique_promoter = with_unique_promoter / with_coordinate_data,
    .groups = "drop"
  )

####################### COMPARATIVE ANALYSIS

# Compare promoter status between CRISPR-covered vs non-covered transcripts

# Get CRISPR-covered transcripts and their promoter status
crispr_covered_transcripts <- data.frame(
  transcript_id = transcripts_with_crispr,
  stringsAsFactors = FALSE
) %>%
  mutate(
    has_unique_promoter = transcript_id %in% unique_promoter_transcripts,
    has_coordinate_data = transcript_id %in% transcript_coords$transcript_id,
    crispr_status = "CRISPR_covered"
  )

# Combine with non-overlapping transcripts
comparison_data <- bind_rows(
  non_overlap_with_promoter_info %>% 
    select(transcript_id, has_unique_promoter, has_coordinate_data) %>%
    mutate(crispr_status = "Non_CRISPR"),
  crispr_covered_transcripts %>% 
    select(transcript_id, has_unique_promoter, has_coordinate_data, crispr_status)
)

# Summary by group
comparison_summary <- comparison_data %>%
  filter(has_coordinate_data) %>%  # Only include transcripts with coordinate data
  group_by(crispr_status) %>%
  summarise(
    total = n(),
    unique_promoter = sum(has_unique_promoter, na.rm = TRUE),
    overlapping_promoter = sum(!has_unique_promoter, na.rm = TRUE),
    prop_unique = unique_promoter / total,
    .groups = "drop"
  )

print("Comparison of promoter status between CRISPR-covered and non-covered transcripts:")
print(comparison_summary)

####################### VISUALIZATION

# Create visualization of promoter status comparison
plot_data <- comparison_data %>%
  filter(has_coordinate_data) %>%
  mutate(
    promoter_status = ifelse(has_unique_promoter, "Unique Promoter", "Overlapping Promoter"),
    crispr_label = ifelse(crispr_status == "CRISPR_covered", "CRISPR Covered", "Non-CRISPR")
  )

write.csv(plot_data, "unique_promoters/transcripts_non_screen_unique_promoter.csv", row.names = FALSE)

# Bar plot showing proportions
ggplot(plot_data, aes(x = crispr_label, fill = promoter_status)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Unique Promoter" = "lightblue", "Overlapping Promoter" = "lightcoral")) +
  labs(
    title = "Promoter Status: CRISPR-Covered vs Non-CRISPR Discordant Transcripts",
    x = "CRISPR Coverage Status",
    y = "Proportion",
    fill = "Promoter Status"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format())

# Count plot
ggplot(plot_data, aes(x = crispr_label, fill = promoter_status)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Unique Promoter" = "lightblue", "Overlapping Promoter" = "lightcoral")) +
  labs(
    title = "Count of Transcripts by Promoter Status and CRISPR Coverage",
    x = "CRISPR Coverage Status",
    y = "Number of Transcripts",
    fill = "Promoter Status"
  ) +
  theme_minimal()
