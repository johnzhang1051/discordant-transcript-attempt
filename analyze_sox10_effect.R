####################### Setting up environment
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)

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

screen_sequence_map <- readr::read_csv(
  file = file.path("depmap-data", "ScreenSequenceMap.csv")
)

####################### FILTER FOR SOX10 GUIDES

# Get all guides targeting SOX10
sox10_guides <- avana_guide_map %>%
  filter(Gene == "SOX10 (6663)") %>%
  select(sgRNA, Gene, GenomeAlignment)

cat("Found", nrow(sox10_guides), "guides targeting SOX10\n")

####################### GET SOX10 GUIDE EFFECTS

# Filter log fold change data for SOX10 guides
sox10_lfc <- avana_logfold_change %>%
  filter(sgRNA %in% sox10_guides$sgRNA) %>%
  pivot_longer(cols = -sgRNA, 
               names_to = "SequenceID", 
               values_to = "log_fold_change")

# Join with cell line metadata
sox10_effects <- sox10_lfc %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID")
  #filter(!is.na(OncotreeLineage))  # Remove any without cancer type info

####################### DEFINE MELANOMA CELL LINES

# Filter to melanoma types (excluding uveal, acral, mucosal)
melanoma_subtypes <- c("Melanoma", "Cutaneous Melanoma")
sox10_effects <- sox10_effects %>%
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

####################### SUMMARIZE BY GUIDE AND CANCER TYPE

# Aggregate effects by guide and cancer lineage
sox10_summary <- sox10_effects %>%
  group_by(sgRNA, OncotreePrimaryDisease) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    median_lfc = median(log_fold_change, na.rm = TRUE),
    sd_lfc = sd(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

####################### PLOT: SOX10 Guide Effects by Cancer Type

# Get the top 25 cancer types with lowest (most negative) mean LFC
top_25_cancers <- sox10_summary %>%
  group_by(OncotreePrimaryDisease) %>%
  summarise(overall_mean_lfc = mean(mean_lfc, na.rm = TRUE)) %>%
  arrange(overall_mean_lfc) %>%
  head(25) %>%
  pull(OncotreePrimaryDisease)

# Filter to only these top 25
sox10_primary_filtered <- sox10_primary_summary %>%
  filter(OncotreePrimaryDisease %in% top_25_cancers)%>%
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
  # If melanoma didn't make top 25, just use gray for all
  color_palette <- setNames(rep("gray70", length(top_25_cancers)), top_25_cancers)
}

ggplot(sox10_primary_filtered, aes(x = reorder(OncotreePrimaryDisease, -mean_lfc), 
                                   y = mean_lfc, 
                                   fill = OncotreePrimaryDisease)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "SOX10 Guide Effects Across Cancer Types (Top 25 Most Negative)",
    subtitle = "Lower values = more cell death (more essential)",
    x = "Cancer Type",
    y = "Mean Log2 Fold Change",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = color_palette)

####################### PLOT 2: Individual Guides - Melanoma vs Others

# Compare melanoma vs non-melanoma for each guide
sox10_melanoma_comparison <- sox10_effects %>%
  mutate(cancer_group = ifelse(is_melanoma, "Melanoma", "Other Cancers")) %>%
  group_by(sgRNA, cancer_group) %>%
  summarise(
    mean_lfc = mean(log_fold_change, na.rm = TRUE),
    n_cell_lines = n(),
    .groups = "drop"
  )

ggplot(sox10_melanoma_comparison, aes(x = sgRNA, y = mean_lfc, fill = cancer_group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  labs(
    title = "SOX10 Guide Effects: Melanoma vs Other Cancers",
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
melanoma_stats <- sox10_effects %>%
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
melanoma_lfc <- sox10_effects %>% filter(is_melanoma) %>% pull(log_fold_change)
other_lfc <- sox10_effects %>% filter(!is_melanoma) %>% pull(log_fold_change)
t_test_result <- t.test(melanoma_lfc, other_lfc)

cat("p-value:", t_test_result$p.value, "\n")
cat("Mean difference:", t_test_result$estimate[1] - t_test_result$estimate[2], "\n")

###################### ANALYZE SOX10 TRANSCRIPTS AND GUIDE EFFECTS

library(biomaRt)
library(GenomicRanges)

# Query biomaRt for SOX10 transcripts
ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get all SOX10 transcripts with their exon information
sox10_transcripts <- getBM(
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
  values = "SOX10",
  mart = ensembl_data
)

cat("\nFound", length(unique(sox10_transcripts$ensembl_transcript_id)), "SOX10 transcripts\n")

####################### MAP GUIDES TO SOX10 TRANSCRIPTS

# Prepare guide data
sox10_guides_map <- sox10_guides %>%
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
sox10_exons <- sox10_transcripts %>%
  mutate(strand_char = case_when(
    strand == 1 ~ "+",
    strand == -1 ~ "-",
    TRUE ~ "*"
  )) %>%
  dplyr::select(-strand)  # Remove conflicting column names

# Create GenomicRanges objects
exon_gr <- makeGRangesFromDataFrame(
  sox10_exons,
  seqnames.field = "chromosome_name",
  start.field = "exon_chrom_start", 
  end.field = "exon_chrom_end",
  strand.field = "strand_char",
  keep.extra.columns = TRUE
  , ignore.strand = TRUE
)

guide_gr <- makeGRangesFromDataFrame(
  sox10_guides_map %>%
    mutate(strand_char = case_when(
      guide_strand == 1 ~ "+",
      guide_strand == -1 ~ "-", 
      TRUE ~ "*"
    )),
  seqnames.field = "chromosome",
  start.field = "guide_coordinate",
  end.field = "guide_coordinate",
  strand.field = "strand_char",
  keep.extra.columns = TRUE
  , ignore.strand = TRUE
)

# Find overlaps between guides and exons
overlap_hits <- findOverlaps(guide_gr, exon_gr)

# Extract overlapping data
sox10_overlaps <- data.frame(
  guide_data = sox10_guides_map[queryHits(overlap_hits), ],
  exon_data = sox10_exons[subjectHits(overlap_hits), ]
)

####################### TRANSCRIPT-LEVEL GUIDE SUMMARY

# Count guides per transcript
sox10_transcript_guide_counts <- sox10_overlaps %>%
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
all_sox10_transcripts <- unique(sox10_transcripts$ensembl_transcript_id)
transcripts_with_guides <- unique(sox10_overlaps$exon_data.ensembl_transcript_id)
transcripts_without_guides <- setdiff(all_sox10_transcripts, transcripts_with_guides)

sox10_transcripts_no_guides <- sox10_transcripts %>%
  filter(ensembl_transcript_id %in% transcripts_without_guides) %>%
  distinct(ensembl_transcript_id, transcript_biotype) %>%
  mutate(
    transcript_id = ensembl_transcript_id,
    n_guides_targeting = 0,
    n_exons_targeted = 0,
    guide_list = ""
  )

# Combine
sox10_all_transcript_summary <- bind_rows(
  sox10_transcript_guide_counts %>% select(-exon_data.ensembl_transcript_id),
  sox10_transcripts_no_guides %>% select(-ensembl_transcript_id)
)

####################### GET GUIDE EFFECTS FOR EACH TRANSCRIPT

# Join overlaps with guide effects
sox10_transcript_effects <- sox10_overlaps %>%
  left_join(
    sox10_lfc %>% rename(guide_sgRNA = sgRNA),
    by = c("guide_data.sgRNA" = "guide_sgRNA")
  ) %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID")

# Define melanoma cell lines
sox10_transcript_effects <- sox10_transcript_effects %>%
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

# Aggregate guide effects by transcript
sox10_transcript_aggregated <- sox10_transcript_effects %>%
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
sox10_transcript_melanoma <- sox10_transcript_aggregated %>%
  filter(is_melanoma) %>%
  select(transcript_id, biotype, n_guides, 
         melanoma_mean_lfc = mean_lfc,
         melanoma_median_lfc = median_lfc,
         melanoma_n_cell_lines = n_cell_lines)

sox10_transcript_other <- sox10_transcript_aggregated %>%
  filter(!is_melanoma) %>%
  select(transcript_id,
         other_mean_lfc = mean_lfc,
         other_median_lfc = median_lfc,
         other_n_cell_lines = n_cell_lines)

# Combine for comparison
sox10_transcript_comparison <- sox10_transcript_melanoma %>%
  left_join(sox10_transcript_other, by = "transcript_id") %>%
  mutate(
    lfc_difference = melanoma_mean_lfc - other_mean_lfc,
    melanoma_specific = melanoma_mean_lfc < -0.5 & other_mean_lfc > -0.3
  )

write.csv(sox10_transcript_comparison, "sox10_transcript_effects_comparison.csv", row.names = FALSE)

####################### TRANSCRIPT VISUALIZATIONS

# Plot 1: Number of guides targeting each transcript
ggplot(sox10_all_transcript_summary, aes(x = reorder(transcript_id, -n_guides_targeting), 
                                         y = n_guides_targeting)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Number of Guides Targeting Each SOX10 Transcript",
    x = "Transcript ID",
    y = "Number of Guides"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Plot 2: Guide effects - Melanoma vs Other cancers for each transcript
sox10_transcript_comparison_long <- sox10_transcript_comparison %>%
  select(transcript_id, melanoma_mean_lfc, other_mean_lfc) %>%
  pivot_longer(cols = c(melanoma_mean_lfc, other_mean_lfc),
               names_to = "cancer_type",
               values_to = "mean_lfc") %>%
  mutate(cancer_type = ifelse(cancer_type == "melanoma_mean_lfc", 
                              "Melanoma", "Other Cancers"))

ggplot(sox10_transcript_comparison_long, 
       aes(x = reorder(transcript_id, -mean_lfc), y = mean_lfc, fill = cancer_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  labs(
    title = "SOX10 Transcript Guide Effects: Melanoma vs Other Cancers",
    subtitle = "Only transcripts with guides shown",
    x = "Transcript ID",
    y = "Mean Log2 Fold Change",
    fill = "Cancer Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("Melanoma" = "red", "Other Cancers" = "gray70"))

# Plot 3: Melanoma specificity - difference in effect
ggplot(sox10_transcript_comparison, 
       aes(x = reorder(transcript_id, -lfc_difference), 
           y = lfc_difference,
           fill = melanoma_specific)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Melanoma Specificity of SOX10 Transcript Targeting",
    subtitle = "Difference = Melanoma LFC - Other Cancers LFC (more negative = more melanoma-specific)",
    x = "Transcript ID",
    y = "LFC Difference (Melanoma - Other)",
    fill = "Melanoma Specific"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  ) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"))

####################### SUMMARY STATISTICS

cat("\n=== SOX10 TRANSCRIPT ANALYSIS SUMMARY ===\n\n")
cat("Total SOX10 transcripts:", length(all_sox10_transcripts), "\n")
cat("Transcripts with guides:", length(transcripts_with_guides), "\n")
cat("Transcripts without guides:", length(transcripts_without_guides), "\n\n")

if (nrow(sox10_transcript_comparison) > 0) {
  melanoma_specific_transcripts <- sox10_transcript_comparison %>%
    filter(melanoma_specific)
  
  cat("MELANOMA-SPECIFIC TRANSCRIPTS (LFC < -0.5 in melanoma, > -0.3 in others):\n")
  if (nrow(melanoma_specific_transcripts) > 0) {
    print(melanoma_specific_transcripts %>% 
            select(transcript_id, n_guides, melanoma_mean_lfc, other_mean_lfc, lfc_difference))
  } else {
    cat("None found with current thresholds.\n")
  }
  
  cat("\n\nALL TRANSCRIPT EFFECTS (sorted by melanoma specificity):\n")
  print(sox10_transcript_comparison %>% 
          arrange(lfc_difference) %>%
          select(transcript_id, biotype, n_guides, melanoma_mean_lfc, other_mean_lfc, lfc_difference))
}

# Export final summary
write.csv(sox10_melanoma_comparison, "sox10_melanoma_comparison.csv", row.names = FALSE)
write.csv(melanoma_stats, "sox10_statistical_summary.csv", row.names = FALSE)

cat("\nAll analysis complete! Files saved.\n")
