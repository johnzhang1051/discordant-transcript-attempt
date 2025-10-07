####################### Setting up environment
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)
library(biomaRt)
library(GenomicRanges)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::setdiff)

####################### LOAD DATA
# Load the guide mapping and effect data
avana_guide_map <- readr::read_csv(
  file=file.path("depmap-data", file="AvanaGuideMap.csv")
)

avana_logfold_change <- readr::read_csv(
  file=file.path("depmap-data", file="AvanaLogfoldChange.csv")
)
colnames(avana_logfold_change)[1] <- "sgRNA"

cell_info <- read.csv("depmap-data/Model.csv")
screen_sequence_map <- readr::read_csv(
  file=file.path("depmap-data", file="ScreenSequenceMap.csv")
)

####################### FILTER FOR EIF3E GUIDES

# Get all guides targeting eif3e
eif3e_guides <- avana_guide_map %>%
  filter(Gene == "EIF3E (3646)") %>%
  select(sgRNA, Gene, GenomeAlignment)

####################### GET EIF3E GUIDE EFFECTS

# Filter log fold change data for eif3e guides
eif3e_lfc <- avana_logfold_change %>%
  filter(sgRNA %in% eif3e_guides$sgRNA) %>%
  pivot_longer(cols = -sgRNA, 
               names_to = "SequenceID", 
               values_to = "log_fold_change")

# Join with cell line metadata
eif3e_effects <- eif3e_lfc %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID")

####################### DEFINE MELANOMA CELL LINES

# Filter to melanoma types (excluding uveal, acral, mucosal)
melanoma_subtypes <- c("Melanoma", "Cutaneous Melanoma")
eif3e_effects <- eif3e_effects %>%
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

####################### SUMMARIZE BY GUIDE AND CANCER TYPE

# Aggregate effects by guide and cancer lineage
eif3e_summary <- eif3e_effects %>%
  group_by(sgRNA, OncotreePrimaryDisease) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    median_lfc = median(log_fold_change, na.rm = TRUE),
    sd_lfc = sd(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

####################### PLOT: EIF3E Guide Effects by Cancer Type

# Get the top 25 cancer types with lowest (most negative) mean LFC
top_25_cancers <- eif3e_summary %>%
  group_by(OncotreePrimaryDisease) %>%
  summarise(overall_mean_lfc = mean(mean_lfc, na.rm = TRUE)) %>%
  arrange(overall_mean_lfc) %>%
  head(25) %>%
  pull(OncotreePrimaryDisease)

# Filter to only these top 25
eif3e_primary_filtered <- eif3e_summary %>%
  filter(OncotreePrimaryDisease %in% top_25_cancers) %>%
  filter(!is.na(OncotreePrimaryDisease), !is.na(mean_lfc))

# Create color palette - check if "Melanoma" is in the filtered data
melanoma_present <- "Melanoma" %in% top_25_cancers
other_cancers <- setdiff(top_25_cancers, "Melanoma")

if (melanoma_present) {
  color_palette <- c(
    "Melanoma" = "red",
    setNames(rep("gray70", length(other_cancers)), other_cancers)
  )
} else {
  color_palette <- setNames(rep("gray70", length(top_25_cancers)), top_25_cancers)
}

ggplot(eif3e_primary_filtered, aes(x = reorder(OncotreePrimaryDisease, -mean_lfc), 
                                   y = mean_lfc, 
                                   fill = OncotreePrimaryDisease)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "EIF3E Guide Effects Across Cancer Types (Top 25 Most Negative)",
    subtitle = "Lower values = more cell death (more essential)",
    x = "Cancer Type",
    y = "Mean Log2 Fold Change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = color_palette)

####################### PLOT 2: Individual Guides - Melanoma vs Others

# Compare melanoma vs non-melanoma for each guide
eif3e_melanoma_comparison <- eif3e_effects %>%
  mutate(cancer_group = ifelse(is_melanoma, "Melanoma", "Other Cancers")) %>%
  group_by(sgRNA, cancer_group) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

ggplot(eif3e_melanoma_comparison, aes(x = sgRNA, y = mean_lfc, fill = cancer_group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  labs(
    title = "EIF3E Guide Effects: Melanoma vs Other Cancers",
    subtitle = "Each bar represents mean effect across cell lines",
    x = "sgRNA",
    y = "Mean Log2 Fold Change",
    fill = "Cancer Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other Cancers" = "gray70"))

####################### STATISTICAL SIGNIFICANCE

# Calculate mean effect in melanoma vs others
melanoma_stats <- eif3e_effects %>%
  mutate(cancer_group = ifelse(is_melanoma, "Melanoma", "Other")) %>%
  group_by(cancer_group) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    median_lfc = median(log_fold_change, na.rm = TRUE),
    sd_lfc = sd(log_fold_change, na.rm = TRUE),
    n = n()
  )

print(melanoma_stats)

# T-test comparing melanoma vs others
melanoma_lfc <- eif3e_effects %>% filter(is_melanoma) %>% pull(log_fold_change)
other_lfc <- eif3e_effects %>% filter(!is_melanoma) %>% pull(log_fold_change)
t_test_result <- t.test(melanoma_lfc, other_lfc)

cat("p-value:", t_test_result$p.value, "\n")
cat("Mean difference:", t_test_result$estimate[1] - t_test_result$estimate[2], "\n")

###################### ANALYZE EIF3E TRANSCRIPTS AND GUIDE EFFECTS

# Query biomaRt for eif3e transcripts
ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get all eif3e transcripts with their exon information
eif3e_transcripts <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "external_gene_name",
    "ensembl_exon_id",
    "exon_chrom_start",
    "exon_chrom_end",
    "strand",
    "chromosome_name",
    "transcript_biotype"
  ),
  filters = "external_gene_name",
  values = "EIF3E",
  mart = ensembl_data
)

cat("\nFound", length(unique(eif3e_transcripts$ensembl_transcript_id)), "EIF3E transcripts\n")

####################### MAP GUIDES TO EIF3E TRANSCRIPTS

# Prepare guide data
eif3e_guides_map <- eif3e_guides %>%
  filter(grepl("^chr[0-9XY]", GenomeAlignment)) %>%
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

# Prepare exon data for GenomicRanges
eif3e_exons <- eif3e_transcripts %>%
  mutate(strand_char = case_when(
    strand == 1 ~ "+",
    strand == -1 ~ "-",
    TRUE ~ "*"
  )) %>%
  dplyr::select(-strand)

# Create GenomicRanges objects
exon_gr <- makeGRangesFromDataFrame(
  eif3e_exons,
  seqnames.field = "chromosome_name",
  start.field = "exon_chrom_start", 
  end.field = "exon_chrom_end",
  strand.field = "strand_char",
  keep.extra.columns = TRUE,
  ignore.strand = TRUE
)

guide_gr <- makeGRangesFromDataFrame(
  eif3e_guides_map %>%
    mutate(strand_char = case_when(
      guide_strand == 1 ~ "+",
      guide_strand == -1 ~ "-", 
      TRUE ~ "*"
    )),
  seqnames.field = "chromosome",
  start.field = "guide_coordinate",
  end.field = "guide_coordinate",
  strand.field = "strand_char",
  keep.extra.columns = TRUE,
  ignore.strand = TRUE
)

# Find overlaps between guides and exons
overlap_hits <- findOverlaps(guide_gr, exon_gr)

# Extract overlapping data
eif3e_overlaps <- data.frame(
  guide_data = eif3e_guides_map[queryHits(overlap_hits), ],
  exon_data = eif3e_exons[subjectHits(overlap_hits), ]
)

####################### TRANSCRIPT-LEVEL GUIDE SUMMARY

# Count guides per transcript
eif3e_transcript_guide_counts <- eif3e_overlaps %>%
  group_by(exon_data.ensembl_transcript_id) %>%
  summarise(
    transcript_id = first(exon_data.ensembl_transcript_id),
    transcript_biotype = first(exon_data.transcript_biotype),
    n_guides_targeting = n_distinct(guide_data.sgRNA),
    n_exons_targeted = n_distinct(exon_data.ensembl_exon_id),
    guide_list = paste(unique(guide_data.sgRNA), collapse = ";"),
    .groups = "drop"
  )

# Get transcripts WITHOUT guides
all_eif3e_transcripts <- unique(eif3e_transcripts$ensembl_transcript_id)
transcripts_with_guides <- unique(eif3e_overlaps$exon_data.ensembl_transcript_id)
transcripts_without_guides <- setdiff(all_eif3e_transcripts, transcripts_with_guides)

eif3e_transcripts_no_guides <- eif3e_transcripts %>%
  filter(ensembl_transcript_id %in% transcripts_without_guides) %>%
  distinct(ensembl_transcript_id, transcript_biotype) %>%
  mutate(
    transcript_id = ensembl_transcript_id,
    n_guides_targeting = 0,
    n_exons_targeted = 0,
    guide_list = ""
  )

# Combine
eif3e_all_transcript_summary <- bind_rows(
  eif3e_transcript_guide_counts %>% select(-exon_data.ensembl_transcript_id),
  eif3e_transcripts_no_guides %>% select(-ensembl_transcript_id)
)

####################### GET GUIDE EFFECTS FOR EACH TRANSCRIPT

# Join overlaps with guide effects
eif3e_transcript_effects <- eif3e_overlaps %>%
  left_join(
    eif3e_lfc %>% rename(guide_sgRNA = sgRNA),
    by = c("guide_data.sgRNA" = "guide_sgRNA")
  ) %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID")

# Define melanoma cell lines
eif3e_transcript_effects <- eif3e_transcript_effects %>%
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

# Aggregate guide effects by transcript
eif3e_transcript_aggregated <- eif3e_transcript_effects %>%
  filter(!is.na(log_fold_change)) %>%
  group_by(exon_data.ensembl_transcript_id, is_melanoma) %>%
  summarise(
    transcript_id = first(exon_data.ensembl_transcript_id),
    biotype = first(exon_data.transcript_biotype),
    n_guides = n_distinct(guide_data.sgRNA),
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    median_lfc = median(log_fold_change, na.rm = TRUE),
    sd_lfc = sd(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

# Separate melanoma vs other cancers
eif3e_transcript_melanoma <- eif3e_transcript_aggregated %>%
  filter(is_melanoma) %>%
  select(transcript_id, biotype, n_guides, 
         melanoma_mean_lfc = mean_lfc,
         melanoma_median_lfc = median_lfc,
         melanoma_n_cell_lines = n_cell_lines)

eif3e_transcript_other <- eif3e_transcript_aggregated %>%
  filter(!is_melanoma) %>%
  select(transcript_id,
         other_mean_lfc = mean_lfc,
         other_median_lfc = median_lfc,
         other_n_cell_lines = n_cell_lines)

# Combine for comparison
eif3e_transcript_comparison <- eif3e_transcript_melanoma %>%
  left_join(eif3e_transcript_other, by = "transcript_id") %>%
  mutate(
    effect_size = melanoma_mean_lfc - other_mean_lfc
  ) %>%
  # Add guide list and cell line info
  left_join(
    eif3e_transcript_guide_counts %>% select(transcript_id, guide_list, n_exons_targeted),
    by = "transcript_id"
  ) %>%
  # Get melanoma cell lines for this transcript
  left_join(
    eif3e_transcript_effects %>%
      filter(is_melanoma) %>%
      group_by(exon_data.ensembl_transcript_id) %>%
      summarise(
        melanoma_cell_lines = paste(unique(ModelID), collapse = ";"),
        .groups = "drop"
      ),
    by = c("transcript_id" = "exon_data.ensembl_transcript_id")
  ) %>%
  # Get non-melanoma cell lines
  left_join(
    eif3e_transcript_effects %>%
      filter(!is_melanoma) %>%
      group_by(exon_data.ensembl_transcript_id) %>%
      summarise(
        non_melanoma_cell_lines = paste(unique(ModelID), collapse = ";"),
        .groups = "drop"
      ),
    by = c("transcript_id" = "exon_data.ensembl_transcript_id")
  ) %>%
  filter(biotype == "protein_coding") %>%
  arrange(effect_size)


write.csv(eif3e_transcript_comparison, "eif3e_analysis/eif3e_transcript_effects_comparison.csv", row.names = FALSE)

####################### COMPARE GUIDES: DISCORDANT vs OTHER (GUIDE-LEVEL ANALYSIS)

# Load discordant transcript list
discordant_list <- read.csv("resubmission_data/discordant_RESUBMISSION.csv")

# Find which EIF3E transcripts are discordant
eif3e_discordant <- eif3e_transcript_comparison %>%
  filter(transcript_id %in% discordant_list$transcript_id) %>%
  pull(transcript_id)

# All other EIF3E transcripts
eif3e_other <- eif3e_transcript_comparison %>%
  filter(!transcript_id %in% discordant_list$transcript_id) %>%
  pull(transcript_id)

cat("\n=== EIF3E TRANSCRIPT CLASSIFICATION ===\n")
cat("Discordant EIF3E transcripts:", length(eif3e_discordant), "\n")
cat("Other EIF3E transcripts:", length(eif3e_other), "\n")

if (length(eif3e_discordant) > 0) {
  cat("\nDiscordant:", paste(eif3e_discordant, collapse = ", "), "\n")
}

# Get guides targeting discordant vs other transcripts
guides_discordant <- eif3e_overlaps %>%
  filter(exon_data.ensembl_transcript_id %in% eif3e_discordant) %>%
  pull(guide_data.sgRNA) %>%
  unique()

guides_other <- eif3e_overlaps %>%
  filter(exon_data.ensembl_transcript_id %in% eif3e_other) %>%
  pull(guide_data.sgRNA) %>%
  unique()

cat("\nGuides targeting discordant transcripts:", length(guides_discordant), "\n")
cat("Guides targeting other transcripts:", length(guides_other), "\n")

####################### GUIDE-LEVEL COMPARISON (NO TRANSCRIPT INFLATION)

# Classify guides at the GUIDE level (not inflated by transcript counts)
eif3e_guide_effects <- eif3e_effects %>%
  mutate(
    transcript_category = case_when(
      sgRNA %in% guides_discordant ~ "Discordant",
      sgRNA %in% guides_other ~ "Other",
      TRUE ~ "Unclassified"
    )
  ) 

# Aggregate by GUIDE, not by (guide × transcript)
guide_level_comparison <- eif3e_guide_effects %>%
  group_by(sgRNA, transcript_category, is_melanoma) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

# Overall comparison across all cell lines
overall_by_guide <- guide_level_comparison %>%
  group_by(transcript_category) %>%
  summarise(
    n_guides = n_distinct(sgRNA),
    mean_lfc = mean(mean_lfc, na.rm = TRUE),
    median_lfc = median(mean_lfc, na.rm = TRUE),
    sd_lfc = sd(mean_lfc, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== GUIDE-LEVEL COMPARISON (ALL CELL LINES) ===\n")
print(overall_by_guide)

# T-test at guide level (all cell lines)
discordant_guides_lfc <- guide_level_comparison %>%
  filter(transcript_category == "Discordant") %>%
  pull(mean_lfc)

other_guides_lfc <- guide_level_comparison %>%
  filter(transcript_category == "Other") %>%
  pull(mean_lfc)

cat("\nNumber of guide measurements (Discordant):", length(discordant_guides_lfc), "\n")
cat("Number of guide measurements (Other):", length(other_guides_lfc), "\n")

####################### ANALYSIS 2: Split by Melanoma vs Non-Melanoma

# Comparison split by cancer type at guide level
cancer_comparison_guides <- guide_level_comparison %>%
  group_by(transcript_category, is_melanoma) %>%
  summarise(
    n_guides = n_distinct(sgRNA),
    mean_lfc = mean(mean_lfc, na.rm = TRUE),
    median_lfc = median(mean_lfc, na.rm = TRUE),
    sd_lfc = sd(mean_lfc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cancer_type = ifelse(is_melanoma, "Melanoma", "Non-Melanoma"))

cat("\n=== GUIDE-LEVEL COMPARISON BY CANCER TYPE ===\n")
print(cancer_comparison_guides)

# T-test for melanoma (guide-level)
discordant_melanoma_guides <- guide_level_comparison %>%
  filter(transcript_category == "Discordant", is_melanoma) %>%
  pull(mean_lfc)

other_melanoma_guides <- guide_level_comparison %>%
  filter(transcript_category == "Other", is_melanoma) %>%
  pull(mean_lfc)

cat("\nMelanoma - Discordant guides:", length(discordant_melanoma_guides), "\n")
cat("Melanoma - Other guides:", length(other_melanoma_guides), "\n")

# Visualization: Guide-level comparison
ggplot(cancer_comparison_guides, 
       aes(x = transcript_category, y = mean_lfc, fill = cancer_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  labs(
    title = "EIF3E Guide Effects: Discordant vs Other Transcripts (Guide-Level)",
    subtitle = "Split by cancer type",
    x = "Transcript Category",
    y = "Mean Log2 Fold Change",
    fill = "Cancer Type"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("Melanoma" = "red", "Non-Melanoma" = "gray70"))

# Export results
write.csv(guide_level_comparison, 
          "eif3e_analysis/discordant_vs_other_guide_level.csv", 
          row.names = FALSE)

write.csv(guide_level_comparison, 
          "eif3e_analysis/discordant_vs_other_guide_level.csv", 
          row.names = FALSE)

