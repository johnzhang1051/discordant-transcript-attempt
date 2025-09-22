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

####################### UNIQUE PROMOTER ANALYSIS WITH PROMOTER COUNTS

# Combine coordinates from all three datasets
all_transcript_coords <- bind_rows(
  Kenny %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Kenny"),
  Laurette %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Laurette"), 
  Louph %>% select(transcript_id, seqnames, start, end, strand, geneId) %>% mutate(source = "Louph")
) %>%
  distinct(transcript_id, .keep_all = TRUE)  # Keep first occurrence if transcript appears in multiple datasets

# Extract transcript coordinates and define promoter regions
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
overlap_df <- as.data.frame(overlaps_promoter) %>%
  filter(queryHits != subjectHits)  # Exclude self-overlaps

# Count how many total promoters are in each promoter region (including self)
promoter_counts <- overlap_df %>%
  group_by(queryHits) %>%
  summarise(n_other_promoters = n(), .groups = "drop") %>%
  mutate(
    transcript_id = transcript_coords$transcript_id[queryHits],
    n_promoters = n_other_promoters + 1  # Add 1 to include the transcript itself
  ) %>%
  select(transcript_id, n_promoters)

# Create comprehensive promoter annotation
promoter_annotation <- transcript_coords %>%
  select(transcript_id) %>%
  left_join(promoter_counts, by = "transcript_id") %>%
  mutate(
    n_promoters = ifelse(is.na(n_promoters), 1, n_promoters),  # Transcripts with no overlaps have 1 promoter
    has_unique_promoter = n_promoters == 1,
    promoter_category = case_when(
      n_promoters == 1 ~ "Unique (1 promoter)",
      n_promoters == 2 ~ "2 promoters",
      n_promoters <= 5 ~ "3-5 promoters",
      n_promoters <= 10 ~ "6-10 promoters",
      TRUE ~ ">10 promoters"
    )
  )

####################### ANALYZE NON-OVERLAPPING DISCORDANT TRANSCRIPTS WITH PROMOTER COUNTS

# Analyze non-overlapping transcripts for promoter information
non_overlap_with_promoter_info <- non_overlap_transcripts %>%
  left_join(promoter_annotation, by = "transcript_id") %>%
  mutate(
    has_coordinate_data = !is.na(n_promoters),
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

# Detailed promoter category breakdown
promoter_category_summary <- non_overlap_with_promoter_info %>%
  filter(has_coordinate_data) %>%
  count(promoter_category) %>%
  mutate(proportion = n / sum(n))

####################### COMPARATIVE ANALYSIS WITH PROMOTER COUNTS

# Compare promoter counts between CRISPR-covered vs non-covered transcripts
crispr_covered_with_promoter <- data.frame(
  transcript_id = transcripts_with_crispr,
  stringsAsFactors = FALSE
) %>%
  left_join(promoter_annotation, by = "transcript_id") %>%
  mutate(
    has_coordinate_data = !is.na(n_promoters),
    crispr_status = "CRISPR_covered"
  )

# Combine datasets for comparison
comparison_data <- bind_rows(
  non_overlap_with_promoter_info %>% 
    select(transcript_id, n_promoters, has_unique_promoter, has_coordinate_data, promoter_category) %>%
    mutate(crispr_status = "Non_CRISPR"),
  crispr_covered_with_promoter %>% 
    select(transcript_id, n_promoters, has_unique_promoter, has_coordinate_data, promoter_category, crispr_status)
)

# Summary by group including promoter counts
comparison_summary_detailed <- comparison_data %>%
  filter(has_coordinate_data) %>%
  group_by(crispr_status) %>%
  summarise(
    total = n(),
    unique_promoter = sum(has_unique_promoter, na.rm = TRUE),
    shared_promoter = sum(!has_unique_promoter, na.rm = TRUE),
    prop_unique = unique_promoter / total,
    mean_promoters = mean(n_promoters, na.rm = TRUE),
    median_promoters = median(n_promoters, na.rm = TRUE),
    max_promoters = max(n_promoters, na.rm = TRUE),
    .groups = "drop"
  )

print("Detailed comparison including promoter overlap counts:")
print(comparison_summary_detailed)

####################### VISUALIZATION WITH PROMOTER COUNTS

# Create visualization of promoter category distribution
plot_data_detailed <- comparison_data %>%
  filter(has_coordinate_data) %>%
  mutate(
    crispr_label = ifelse(crispr_status == "CRISPR_covered", "CRISPR Covered", "Non-CRISPR"),
    promoter_category = factor(promoter_category, 
                               levels = c("Unique (0 overlaps)", "1 overlap", "2-5 overlaps", "6-10 overlaps", ">10 overlaps"))
  )

# Stacked bar plot showing promoter category distribution
ggplot(plot_data_detailed, aes(x = crispr_label, fill = promoter_category)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  labs(
    title = "Promoter Overlap Categories: CRISPR-Covered vs Non-CRISPR Discordant Transcripts",
    x = "CRISPR Coverage Status",
    y = "Proportion",
    fill = "Promoter Category"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format())

# Count plot by category
ggplot(plot_data_detailed, aes(x = crispr_label, fill = promoter_category)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  labs(
    title = "Count of Transcripts by Promoter Overlap Category and CRISPR Coverage",
    x = "CRISPR Coverage Status",
    y = "Number of Transcripts",
    fill = "Promoter Category"
  ) +
  theme_minimal()

# Box plot of number of promoters per transcript
plot_data_for_boxplot <- comparison_data %>%
  filter(has_coordinate_data) %>%
  mutate(crispr_label = ifelse(crispr_status == "CRISPR_covered", "CRISPR Covered", "Non-CRISPR"))

ggplot(plot_data_for_boxplot, aes(x = crispr_label, y = n_promoters, fill = crispr_label)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("CRISPR Covered" = "lightblue", "Non-CRISPR" = "lightcoral")) +
  labs(
    title = "Number of Promoters per Transcript",
    x = "CRISPR Coverage Status",
    y = "Number of Promoters",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Save detailed results
write.csv(comparison_data, "unique_promoters/transcripts_promoter_detailed_analysis.csv", row.names = FALSE)
write.csv(promoter_annotation, "unique_promoters/all_transcripts_promoter_data.csv", row.names = FALSE)
