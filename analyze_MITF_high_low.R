####################### Setting up environment
library(depmap)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)

# Fix select() conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

####################### Load Required Data

# Load transcript lists
discordant_transcripts <- readr::read_csv(
  file = file.path("correlation_results", "discordant_transcripts.csv")
)
colnames(discordant_transcripts)[1] <- "transcript_id"

# Load MITF classifications
mitf_classifications <- readr::read_csv(
  file = "mitf_high_low/mitf_expression_classifications_full.csv"
)

# Load overlaps data (transcript-guide mappings from GenomicRanges analysis)
transcript_list_name <- "discordant_transcripts"
overlaps <- readr::read_csv(
  file = paste0("guide_effect/", transcript_list_name, "_exon_guide_overlaps.csv")
)

####################### Load DepMap Data (same as original analysis)

# Load raw read counts for CRISPR screen
avana_raw_readcounts <- readr::read_csv(
  file = file.path("depmap-data", "AvanaRawReadcounts.csv")
)
colnames(avana_raw_readcounts)[1] <- "sg_rna"

screen_sequence_map <- readr::read_csv(
  file = file.path("depmap-data", "ScreenSequenceMap.csv")
)

####################### Process Cell Line Data (Enhanced with MITF)

# Get cell line info and add MITF classification
cell_info <- depmap_metadata()

# Filter to melanoma cells and join with MITF classification
specific_subtypes <- c("Melanoma", "Melanoma, amelanotic")
melanoma_cells <- cell_info %>%
  filter(subtype_disease %in% specific_subtypes) %>%
  left_join(mitf_classifications, by = "depmap_id") %>%
  filter(!is.na(mitf_binary))  # Only keep cells with MITF classification

# Create sequence mapping with MITF status
melanoma_cells_sequences <- melanoma_cells %>% 
  left_join(screen_sequence_map, by = c("depmap_id" = "ModelID")) %>%
  select(SequenceID, subtype_disease.x, depmap_id, mitf_binary, mitf_expression)

####################### Calculate Guide Effects (Same as Original)

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

# Filter to melanoma cells with MITF data
melanoma_effects <- avana_effects %>% 
  inner_join(melanoma_cells_sequences, by = "SequenceID")

# Aggregate guide effects (same as original)
melanoma_effects_aggregated <- melanoma_effects %>%
  group_by(sg_rna) %>%
  summarise(
    mean_read_count = mean(read_count, na.rm = TRUE),
    mean_log_count = mean(log_count, na.rm = TRUE),
    mean_guide_effect = mean(guide_effect, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

####################### Map to Transcripts (Same as Original)

transcript_guide_effects <- overlaps %>% 
  left_join(melanoma_effects_aggregated, by = c("guide_data.sgRNA" = "sg_rna"))

# Average guide effects by transcript (same as original)
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

####################### NOW Split by MITF Status

# Calculate guide effects separately for High and Low MITF
melanoma_effects_by_mitf <- melanoma_effects %>%
  group_by(sg_rna, mitf_binary) %>%
  summarise(
    mean_guide_effect = mean(guide_effect, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

# Map to transcripts by MITF status
transcript_guide_effects_by_mitf <- overlaps %>% 
  left_join(melanoma_effects_by_mitf, by = c("guide_data.sgRNA" = "sg_rna")) %>%
  filter(!is.na(mean_guide_effect))

# Aggregate to transcript level by MITF status
transcript_effects_by_mitf <- transcript_guide_effects_by_mitf %>%
  group_by(exon_data.ensembl_transcript_id, mitf_binary) %>%
  summarise(
    gene_name = dplyr::first(exon_data.external_gene_name),
    n_guides = n(),
    mean_guide_effect = mean(mean_guide_effect, na.rm = TRUE),
    median_guide_effect = median(mean_guide_effect, na.rm = TRUE),
    sd_guide_effect = sd(mean_guide_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_guide_effect))

####################### Create Comparison Table

# Pivot to compare High vs Low MITF effects
comparison_table <- transcript_effects_by_mitf %>%
  select(transcript_id = exon_data.ensembl_transcript_id, 
         gene_name, mitf_binary, mean_guide_effect, n_guides) %>%
  pivot_wider(names_from = mitf_binary, 
              values_from = c(mean_guide_effect, n_guides),
              names_sep = "_") %>%
  # Calculate difference (High MITF - Low MITF)
  mutate(
    effect_difference = mean_guide_effect_High - mean_guide_effect_Low,
    # Classify the difference
    differential_effect = case_when(
      is.na(effect_difference) ~ "Missing data",
      effect_difference < -0.3 ~ "More essential in High MITF",
      effect_difference > 0.3 ~ "More essential in Low MITF", 
      TRUE ~ "Similar effect"
    )
  ) %>%
  arrange(effect_difference)

####################### Analysis and Visualization

# Statistical test if we have paired data
valid_pairs <- comparison_table %>%
  filter(!is.na(mean_guide_effect_High) & !is.na(mean_guide_effect_Low))

if(nrow(valid_pairs) > 10) {
  t_test_result <- t.test(valid_pairs$mean_guide_effect_High, 
                          valid_pairs$mean_guide_effect_Low, 
                          paired = TRUE)
  
  cat("Mean difference (High - Low):", round(t_test_result$estimate, 4), "\n")
  cat("P-value:", format(t_test_result$p.value, scientific = TRUE), "\n")

####################### Visualizations

# 1. Overall distribution (same as original)
ggplot(transcript_aggregated_effects, aes(x = "All Transcripts", y = mean_guide_effect)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(
    title = paste("Guide Effects for", str_to_title(gsub("_", " ", transcript_list_name))),
    x = "",
    y = "Mean Guide Effect",
    subtitle = paste("n =", nrow(transcript_aggregated_effects), "transcripts")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)

# 2. Comparison by MITF status
if(nrow(transcript_effects_by_mitf) > 0) {
  ggplot(transcript_effects_by_mitf, aes(x = mitf_binary, y = mean_guide_effect, fill = mitf_binary)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = c("High" = "darkred", "Low" = "steelblue")) +
    labs(
      title = "Guide Effects Split by MITF Expression",
      x = "MITF Expression Level",
      y = "Mean Guide Effect",
      fill = "MITF Status"
    ) +
    theme_minimal() +
    guides(fill = "none")
}

####################### Export Results

# Create output directory
if (!dir.exists("guide_effect_mitf_split")) {
  dir.create("guide_effect_mitf_split")
}

# Export original overall results
overall_output_file <- paste0("guide_effect_mitf_split/", transcript_list_name, "_guide_effects_overall.csv")
write.csv(transcript_aggregated_effects, overall_output_file, row.names = FALSE)

# Export MITF comparison
if(nrow(comparison_table) > 0) {
  comparison_output_file <- paste0("guide_effect_mitf_split/", transcript_list_name, "_guide_effects_mitf_comparison.csv")
  write.csv(comparison_table, comparison_output_file, row.names = FALSE)
}

# Export detailed by MITF
detailed_output_file <- paste0("guide_effect_mitf_split/", transcript_list_name, "_guide_effects_by_mitf_detailed.csv")
write.csv(transcript_effects_by_mitf, detailed_output_file, row.names = FALSE)


####################### Correlation Analysis: MITF Expression vs Guide Effects
# Calculate guide effects for each cell line individually (not aggregated)
melanoma_effects_by_cell_line <- melanoma_effects %>%
  group_by(sg_rna, depmap_id, mitf_expression) %>%
  summarise(
    mean_guide_effect = mean(guide_effect, na.rm = TRUE),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_guide_effect))

# Map to transcripts for correlation analysis
transcript_guide_effects_by_cell <- overlaps %>% 
  left_join(melanoma_effects_by_cell_line, by = c("guide_data.sgRNA" = "sg_rna")) %>%
  filter(!is.na(mean_guide_effect))

# Calculate correlations for each transcript
transcript_mitf_correlations <- transcript_guide_effects_by_cell %>%
  group_by(exon_data.ensembl_transcript_id) %>%
  filter(n() >= 5) %>%  # Need at least 5 cell lines for meaningful correlation
  summarise(
    gene_name = dplyr::first(exon_data.external_gene_name),
    n_cell_lines = n(),
    n_guides = length(unique(guide_data.sgRNA)),
    pearson_corr = cor(mitf_expression, mean_guide_effect, method = "pearson", use = "complete.obs"),
    spearman_corr = cor(mitf_expression, mean_guide_effect, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) %>%
  filter(!is.na(pearson_corr)) %>%
  # Calculate p-values for correlations
  rowwise() %>%
  mutate(
    # Simple correlation p-value calculation
    pearson_pvalue = if(!is.na(pearson_corr) & n_cell_lines > 3) {
      t_stat <- pearson_corr * sqrt((n_cell_lines - 2) / (1 - pearson_corr^2))
      2 * pt(abs(t_stat), df = n_cell_lines - 2, lower.tail = FALSE)
    } else NA
  ) %>%
  ungroup() %>%
  # Add FDR correction
  mutate(pearson_fdr = p.adjust(pearson_pvalue, method = "BH")) %>%
  arrange(pearson_corr)


