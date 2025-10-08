####################### Setting up environment
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)
library(GenomicRanges)
library(biomaRt)

# Fix select() conflicts
conflict_prefer("select", "dplyr")
conflicts_prefer(GenomicRanges::setdiff)
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::rename)

####################### Load master table with crispr_screened column
master_table <- read_csv("annotated_table/annotated_transcripts.csv")

# Clean transcript IDs for all transcripts (not just non-screened)
master_table$transcript_id_clean <- sub("\\..*", "", sub("^[.]*", "", master_table$transcript_id))

####################### Load existing biomaRt data from saved CSVs
# Load your previously saved biomaRt exon location data
correlated_exon_locations <- read_csv("guide_effect/mitf_high/correlated_MITF_HIGH_exon_locations.csv")
discordant_exon_locations <- read_csv("guide_effect/mitf_high/discordant_MITF_HIGH_exon_locations.csv")

# Combine them
exon_locations <- bind_rows(
  correlated_exon_locations %>% mutate(source_list = "correlated"),
  discordant_exon_locations %>% mutate(source_list = "discordant")
)

####################### 1. COORDINATE OVERLAP ANALYSIS (exon-level overlaps)
cat("1. Analyzing coordinate overlaps at exon level...\n")

# Filter exon_locations to all transcripts in master table
exon_data <- exon_locations %>%
  filter(ensembl_transcript_id %in% master_table$transcript_id_clean) %>%
  filter(!is.na(chromosome_name), !is.na(exon_chrom_start), !is.na(exon_chrom_end))

# Create GenomicRanges object for exons
exon_gr <- GRanges(
  seqnames = exon_data$chromosome_name,
  ranges = IRanges(start = exon_data$exon_chrom_start, 
                   end = exon_data$exon_chrom_end),
  transcript_id = exon_data$ensembl_transcript_id,
  exon_id = exon_data$ensembl_exon_id,
  gene_name = exon_data$external_gene_name
)

# Find overlaps between exons
exon_overlaps <- findOverlaps(exon_gr, exon_gr)
exon_overlap_df <- as.data.frame(exon_overlaps) %>%
  filter(queryHits != subjectHits) %>%  # Exclude self-overlaps
  mutate(
    query_transcript = exon_gr$transcript_id[queryHits],
    subject_transcript = exon_gr$transcript_id[subjectHits]
  ) %>%
  # Only count overlaps between different transcripts
  filter(query_transcript != subject_transcript)

# Count how many other transcripts each transcript overlaps with (at exon level)
transcript_overlap_counts <- exon_overlap_df %>%
  group_by(query_transcript) %>%
  summarise(
    n_overlapping_transcripts = n_distinct(subject_transcript),
    n_overlapping_exons = n(),
    .groups = 'drop'
  ) %>%
  rename(transcript_id = query_transcript)

# Create transcript-level summary with overlap info
transcript_coords_with_overlaps <- exon_data %>%
  group_by(ensembl_transcript_id) %>%
  summarise(
    external_gene_name = first(external_gene_name),
    chromosome_name = first(chromosome_name),
    transcript_start = min(exon_chrom_start, na.rm = TRUE),
    transcript_end = max(exon_chrom_end, na.rm = TRUE),
    transcript_length = transcript_end - transcript_start + 1,
    n_exons = n_distinct(ensembl_exon_id),
    total_exon_length = sum(exon_chrom_end - exon_chrom_start + 1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Join with overlap counts
  left_join(transcript_overlap_counts, by = c("ensembl_transcript_id" = "transcript_id")) %>%
  mutate(
    n_overlapping_transcripts = ifelse(is.na(n_overlapping_transcripts), 0, n_overlapping_transcripts),
    n_overlapping_exons = ifelse(is.na(n_overlapping_exons), 0, n_overlapping_exons),
    has_coordinate_overlaps = n_overlapping_transcripts > 0,
    overlap_category = case_when(
      n_overlapping_transcripts == 0 ~ "No overlaps",
      n_overlapping_transcripts == 1 ~ "Overlaps with 1 transcript",
      n_overlapping_transcripts <= 3 ~ "Overlaps with 2-3 transcripts",
      n_overlapping_transcripts <= 10 ~ "Overlaps with 4-10 transcripts",
      TRUE ~ "Overlaps with >10 transcripts"
    )
  )

####################### GC CONTENT ANALYSIS
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

# Filter to transcripts that have coordinate data
transcripts_with_coords <- transcript_coords_with_overlaps %>%
  filter(!is.na(chromosome_name), !is.na(transcript_start), !is.na(transcript_end))

if (nrow(transcripts_with_coords) > 0) {
  # Create GRanges for sequence extraction
  gr_for_gc <- GRanges(
    seqnames = paste0("chr", transcripts_with_coords$chromosome_name),
    ranges = IRanges(start = transcripts_with_coords$transcript_start, 
                     end = transcripts_with_coords$transcript_end),
    transcript_id = transcripts_with_coords$ensembl_transcript_id
  )
  # Extract sequences for each transcript
  sequences <- tryCatch({
    getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_for_gc)
  }, error = function(e) {
    cat("Error extracting sequences:", e$message, "\n")
    cat("Trying without 'chr' prefix...\n")
    # Try without chr prefix
    gr_for_gc_no_chr <- GRanges(
      seqnames = transcripts_with_coords$chromosome_name,
      ranges = IRanges(start = transcripts_with_coords$transcript_start, 
                       end = transcripts_with_coords$transcript_end),
      transcript_id = transcripts_with_coords$ensembl_transcript_id
    )
    getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_for_gc_no_chr)
  })
  
  if (!is.null(sequences)) {
    # Calculate GC content
    gc_content_values <- letterFrequency(sequences, "GC", as.prob = TRUE)
    
    # Create GC analysis dataframe
    gc_analysis <- data.frame(
      ensembl_transcript_id = transcripts_with_coords$ensembl_transcript_id,
      gene_gc_content = as.numeric(gc_content_values),
      stringsAsFactors = FALSE
    )
  } else {
    gc_analysis <- data.frame(
      ensembl_transcript_id = character(0),
      gene_gc_content = numeric(0),
      stringsAsFactors = FALSE
    )
    cat("Could not extract sequences for GC analysis\n")
  }
} else {
  gc_analysis <- data.frame(
    ensembl_transcript_id = character(0),
    gene_gc_content = numeric(0),
    stringsAsFactors = FALSE
  )
  cat("No transcripts with coordinate data for GC analysis\n")
}

####################### COMBINE ALL ANALYSES
# Merge all analyses - now comparing by screening status instead of correlated/discordant
complete_analysis <- master_table %>%
  # Join with coordinate/overlap analysis
  left_join(transcript_coords_with_overlaps, by = c("transcript_id_clean" = "ensembl_transcript_id")) %>%
  # Join with GC analysis
  left_join(gc_analysis, by = c("transcript_id_clean" = "ensembl_transcript_id")) %>%
  # Create screening status labels
  mutate(
    screening_status = ifelse(crispr_screened, "CRISPR Screened", "Non-CRISPR Screened")
  )

write.csv(complete_analysis, "non_screened_analyses/comparison_analysis_results.csv", row.names = FALSE)

####################### Setting up environment
library(tidyverse)
library(dplyr)
library(ggplot2)
library(conflicted)
library(GenomicRanges)
library(biomaRt)
library(psych)

# Fix select() conflicts
conflict_prefer("select", "dplyr")
conflicts_prefer(GenomicRanges::setdiff)
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::rename)

####################### Load master table with crispr_screened column
master_table <- read_csv("annotated_table/annotated_transcripts.csv")

# Clean transcript IDs for all transcripts (not just non-screened)
master_table$transcript_id_clean <- sub("\\..*", "", sub("^[.]*", "", master_table$transcript_id))

####################### Load existing biomaRt data from saved CSVs
# Load your previously saved biomaRt exon location data
correlated_exon_locations <- read_csv("guide_effect/mitf_high/correlated_MITF_HIGH_exon_locations.csv")
discordant_exon_locations <- read_csv("guide_effect/mitf_high/discordant_MITF_HIGH_exon_locations.csv")

# Combine them
exon_locations <- bind_rows(
  correlated_exon_locations %>% mutate(source_list = "correlated"),
  discordant_exon_locations %>% mutate(source_list = "discordant")
)

####################### 1. COORDINATE OVERLAP ANALYSIS (exon-level overlaps)
cat("1. Analyzing coordinate overlaps at exon level...\n")

# Filter exon_locations to all transcripts in master table
exon_data <- exon_locations %>%
  filter(ensembl_transcript_id %in% master_table$transcript_id_clean) %>%
  filter(!is.na(chromosome_name), !is.na(exon_chrom_start), !is.na(exon_chrom_end))

# Create GenomicRanges object for exons
exon_gr <- GRanges(
  seqnames = exon_data$chromosome_name,
  ranges = IRanges(start = exon_data$exon_chrom_start, 
                   end = exon_data$exon_chrom_end),
  transcript_id = exon_data$ensembl_transcript_id,
  exon_id = exon_data$ensembl_exon_id,
  gene_name = exon_data$external_gene_name
)

# Find overlaps between exons
exon_overlaps <- findOverlaps(exon_gr, exon_gr)
exon_overlap_df <- as.data.frame(exon_overlaps) %>%
  filter(queryHits != subjectHits) %>%  # Exclude self-overlaps
  mutate(
    query_transcript = exon_gr$transcript_id[queryHits],
    subject_transcript = exon_gr$transcript_id[subjectHits]
  ) %>%
  # Only count overlaps between different transcripts
  filter(query_transcript != subject_transcript)

# Count how many other transcripts each transcript overlaps with (at exon level)
transcript_overlap_counts <- exon_overlap_df %>%
  group_by(query_transcript) %>%
  summarise(
    n_overlapping_transcripts = n_distinct(subject_transcript),
    n_overlapping_exons = n(),
    .groups = 'drop'
  ) %>%
  rename(transcript_id = query_transcript)

# Create transcript-level summary with overlap info
transcript_coords_with_overlaps <- exon_data %>%
  group_by(ensembl_transcript_id) %>%
  summarise(
    external_gene_name = first(external_gene_name),
    chromosome_name = first(chromosome_name),
    transcript_start = min(exon_chrom_start, na.rm = TRUE),
    transcript_end = max(exon_chrom_end, na.rm = TRUE),
    transcript_length = transcript_end - transcript_start + 1,
    n_exons = n_distinct(ensembl_exon_id),
    total_exon_length = sum(exon_chrom_end - exon_chrom_start + 1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Join with overlap counts
  left_join(transcript_overlap_counts, by = c("ensembl_transcript_id" = "transcript_id")) %>%
  mutate(
    n_overlapping_transcripts = ifelse(is.na(n_overlapping_transcripts), 0, n_overlapping_transcripts),
    n_overlapping_exons = ifelse(is.na(n_overlapping_exons), 0, n_overlapping_exons),
    has_coordinate_overlaps = n_overlapping_transcripts > 0,
    overlap_category = case_when(
      n_overlapping_transcripts == 0 ~ "No overlaps",
      n_overlapping_transcripts == 1 ~ "Overlaps with 1 transcript",
      n_overlapping_transcripts <= 3 ~ "Overlaps with 2-3 transcripts",
      n_overlapping_transcripts <= 10 ~ "Overlaps with 4-10 transcripts",
      TRUE ~ "Overlaps with >10 transcripts"
    )
  )

####################### GC CONTENT ANALYSIS
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

# Filter to transcripts that have coordinate data
transcripts_with_coords <- transcript_coords_with_overlaps %>%
  filter(!is.na(chromosome_name), !is.na(transcript_start), !is.na(transcript_end))

if (nrow(transcripts_with_coords) > 0) {
  # Create GRanges for sequence extraction
  gr_for_gc <- GRanges(
    seqnames = paste0("chr", transcripts_with_coords$chromosome_name),
    ranges = IRanges(start = transcripts_with_coords$transcript_start, 
                     end = transcripts_with_coords$transcript_end),
    transcript_id = transcripts_with_coords$ensembl_transcript_id
  )
  # Extract sequences for each transcript
  sequences <- tryCatch({
    getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_for_gc)
  }, error = function(e) {
    cat("Error extracting sequences:", e$message, "\n")
    cat("Trying without 'chr' prefix...\n")
    # Try without chr prefix
    gr_for_gc_no_chr <- GRanges(
      seqnames = transcripts_with_coords$chromosome_name,
      ranges = IRanges(start = transcripts_with_coords$transcript_start, 
                       end = transcripts_with_coords$transcript_end),
      transcript_id = transcripts_with_coords$ensembl_transcript_id
    )
    getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_for_gc_no_chr)
  })
  
  if (!is.null(sequences)) {
    # Calculate GC content
    gc_content_values <- letterFrequency(sequences, "GC", as.prob = TRUE)
    
    # Create GC analysis dataframe
    gc_analysis <- data.frame(
      ensembl_transcript_id = transcripts_with_coords$ensembl_transcript_id,
      gene_gc_content = as.numeric(gc_content_values),
      stringsAsFactors = FALSE
    )
  } else {
    gc_analysis <- data.frame(
      ensembl_transcript_id = character(0),
      gene_gc_content = numeric(0),
      stringsAsFactors = FALSE
    )
    cat("Could not extract sequences for GC analysis\n")
  }
} else {
  gc_analysis <- data.frame(
    ensembl_transcript_id = character(0),
    gene_gc_content = numeric(0),
    stringsAsFactors = FALSE
  )
  cat("No transcripts with coordinate data for GC analysis\n")
}

####################### COMBINE ALL ANALYSES
cat("Combining all technical analyses...\n")

# Merge all analyses - now comparing by screening status instead of correlated/discordant
complete_analysis <- master_table %>%
  # Join with coordinate/overlap analysis
  left_join(transcript_coords_with_overlaps, by = c("transcript_id_clean" = "ensembl_transcript_id")) %>%
  # Join with GC analysis
  left_join(gc_analysis, by = c("transcript_id_clean" = "ensembl_transcript_id")) %>%
  # Create screening status labels
  mutate(
    screening_status = ifelse(crispr_screened, "CRISPR Screened", "Non-CRISPR Screened")
  )

####################### COMPARISONS

# OVERLAP COMPARISON - Table of percentages
cat("=== OVERLAP COMPARISON ===\n")
overlap_comparison <- complete_analysis %>%
  filter(!is.na(overlap_category)) %>%
  count(screening_status, overlap_category) %>%
  group_by(screening_status) %>%
  mutate(
    total = sum(n),
    percentage = round(n / total * 100, 1)
  ) %>%
  select(screening_status, overlap_category, n, percentage) %>%
  arrange(screening_status, overlap_category)

print(overlap_comparison)

# Create a cleaner table view
overlap_table <- overlap_comparison %>%
  select(screening_status, overlap_category, percentage) %>%
  pivot_wider(names_from = screening_status, values_from = percentage, values_fill = 0)

cat("\nOverlap Percentage Table:\n")
print(overlap_table)

# Visual comparison - stacked bar chart for overlap categories
overlap_plot <- complete_analysis %>%
  filter(!is.na(overlap_category)) %>%
  count(screening_status, overlap_category) %>%
  group_by(screening_status) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ggplot(aes(x = screening_status, y = percentage, fill = overlap_category)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(percentage > 5, paste0(round(percentage, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  labs(
    title = "Overlap Pattern Distribution: Screened vs Non-Screened Transcripts",
    x = "Screening Status",
    y = "Percentage of Transcripts",
    fill = "Overlap Category"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

print(overlap_plot)

# GC CONTENT COMPARISON - Side-by-side boxplots + summary stats
if (nrow(gc_analysis) > 0 && sum(!is.na(complete_analysis$gene_gc_content)) > 0) {
  
  cat("\n=== GC CONTENT COMPARISON ===\n")
  gc_summary <- complete_analysis %>%
    filter(!is.na(gene_gc_content)) %>%
    group_by(screening_status) %>%
    summarise(
      count = n(),
      mean_gc = round(mean(gene_gc_content), 3),
      median_gc = round(median(gene_gc_content), 3),
      sd_gc = round(sd(gene_gc_content), 3),
      extreme_low_pct = round(mean(gene_gc_content < 0.20) * 100, 1),
      extreme_high_pct = round(mean(gene_gc_content > 0.80) * 100, 1),
      .groups = 'drop'
    )
  
  print(gc_summary)
  
  # Boxplot comparison
  ggplot(complete_analysis %>% filter(!is.na(gene_gc_content)), 
         aes(x = screening_status, y = gene_gc_content, fill = screening_status)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    geom_hline(yintercept = c(0.20, 0.80), linetype = "dashed", color = "red") +
    labs(
      title = "GC Content Comparison: Screened vs Non-Screened",
      x = "Screening Status",
      y = "GC Content (proportion)",
      fill = "Screening Status"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
} else {
  cat("No GC content data available for comparison\n")
}

# TRANSCRIPT LENGTH COMPARISON - Side-by-side boxplots + summary stats
cat("\n=== TRANSCRIPT LENGTH COMPARISON ===\n")
length_summary <- complete_analysis %>%
  filter(!is.na(transcript_length)) %>%
  group_by(screening_status) %>%
  summarise(
    count = n(),
    mean_length = round(mean(transcript_length)),
    median_length = round(median(transcript_length)),
    sd_length = round(sd(transcript_length)),
    min_length = min(transcript_length),
    max_length = max(transcript_length),
    .groups = 'drop'
  )

print(length_summary)

# Boxplot comparison for transcript length
ggplot(complete_analysis %>% filter(!is.na(transcript_length)), 
       aes(x = screening_status, y = log10(transcript_length), fill = screening_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    title = "Transcript Length Comparison: Screened vs Non-Screened",
    x = "Screening Status",
    y = "Log10(Transcript Length)",
    fill = "Screening Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

####################### MITF CORRELATION ANALYSIS
# Analyzing MITF correlation differences between screened and non-screened transcripts

# Use correlation data already in complete_analysis (from master table)
# Filter for transcripts with correlation data
mitf_comparison <- complete_analysis %>%
  filter(!is.na(r_transcript) | !is.na(r_gene))

# Summary statistics by screening status - transcript-level correlations
transcript_mitf_summary <- mitf_comparison %>%
  filter(!is.na(r_transcript)) %>%
  group_by(screening_status) %>%
  summarise(
    count = n(),
    # Transcript-level Pearson correlation stats
    mean_r_transcript = round(mean(r_transcript, na.rm = TRUE), 4),
    median_r_transcript = round(median(r_transcript, na.rm = TRUE), 4),
    sd_r_transcript = round(sd(r_transcript, na.rm = TRUE), 4),
    .groups = 'drop'
  )

# Summary statistics by screening status - gene-level correlations
gene_mitf_summary <- mitf_comparison %>%
  filter(!is.na(r_gene)) %>%
  group_by(screening_status) %>%
  summarise(
    count = n(),
    # Gene-level correlation stats
    mean_gene_corr = round(mean(r_gene, na.rm = TRUE), 4),
    median_gene_corr = round(median(r_gene, na.rm = TRUE), 4),
    sd_gene_corr = round(sd(r_gene, na.rm = TRUE), 4),
    .groups = 'drop'
  )
write.csv(transcript_mitf_summary, "non_screened_analyses/transcript_mitf_correlation_summary.csv", row.names = FALSE)
write.csv(gene_mitf_summary, "non_screened_analyses/gene_mitf_correlation_summary.csv", row.names = FALSE)

####################### UPDATE MASTER TABLE

# add GC % content, Transcript Length, Transcript Coordinate overlaps for each transcript
analytical_data <- complete_analysis %>%
  select(
    transcript_id_clean,
    # Coordinate data
    chromosome_name, transcript_start, transcript_end, transcript_length,
    n_exons, total_exon_length,
    # Overlap data
    n_overlapping_transcripts, n_overlapping_exons, 
    has_coordinate_overlaps, overlap_category,
    # GC content data
    gene_gc_content
  ) %>%
  rename(
    transcript_id = transcript_id_clean,
    gc_content_proportion = gene_gc_content,
    transcript_chr = chromosome_name,
  )

# Update the master table by joining with analytical data
updated_master_table <- master_table %>%
  left_join(analytical_data, by = c("transcript_id_clean" = "transcript_id")) %>%
  # Organize columns logically
  select(
    # Original master table columns first
    transcript_id, transcript_id_clean, transcript_classification, crispr_screened,
    # New GC content information
    gc_content_proportion,
    # Any remaining original columns
    everything(),
    -transcript_id_clean
  ) %>%
  # Sort by classification and transcript_id
  arrange(desc(transcript_classification), desc(crispr_screened), effect_size, guide_effect_fdr)

# Save the updated master table
write.csv(updated_master_table, "annotated_table/annotated_transcripts_updated.csv", row.names = FALSE)


####################### VISUALIZATIONS

# Overlap comparision between screened and non-screened
overlap_comparison <- complete_analysis %>%
  filter(!is.na(overlap_category)) %>%
  count(screening_status, overlap_category) %>%
  group_by(screening_status) %>%
  mutate(
    total = sum(n),
    percentage = n / total
  ) %>%
  select(screening_status, overlap_category, n, percentage) %>%
  arrange(screening_status, overlap_category)

overlap_table <- overlap_comparison %>%
  select(screening_status, overlap_category, percentage) %>%
  pivot_wider(names_from = screening_status, values_from = percentage, values_fill = 0)

View(overlap_table)
write.csv(overlap_table, "non_screened_analyses/overlap_compare_table.csv", row.names = FALSE)

overlap_plot <- complete_analysis %>%
  filter(!is.na(overlap_category)) %>%
  count(screening_status, overlap_category) %>%
  group_by(screening_status) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ggplot(aes(x = screening_status, y = percentage, fill = overlap_category)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(percentage > 5, paste0(round(percentage, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  labs(
    title = "Overlap Distribution: Screened vs Non-Screened Transcripts",
    x = "Screening Status",
    y = "Percentage of Transcripts",
    fill = "Overlap Category"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

print(overlap_plot)

# Comparing GC % content
gc_summary <- complete_analysis %>%
  filter(!is.na(gene_gc_content)) %>%
  group_by(screening_status) %>%
  summarise(
    count = n(),
    mean_gc = round(mean(gene_gc_content), 3),
    median_gc = round(median(gene_gc_content), 3),
    sd_gc = round(sd(gene_gc_content), 3),
    extreme_low_pct = round(mean(gene_gc_content < 0.20) * 100, 1),
    extreme_high_pct = round(mean(gene_gc_content > 0.80) * 100, 1),
    .groups = 'drop'
  )

View(gc_summary)
write.csv(gc_summary, "non_screened_analyses/gc_summary.csv", row.names = FALSE)


ggplot(complete_analysis %>% filter(!is.na(gene_gc_content)), 
       aes(x = screening_status, y = gene_gc_content, fill = screening_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = c(0.20, 0.80), linetype = "dashed", color = "red") +
  labs(
    title = "GC Content Comparison: Screened vs Non-Screened",
    x = "Screening Status",
    y = "GC Content (proportion)",
    fill = "Screening Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Transcript length comparison
length_summary <- complete_analysis %>%
  filter(!is.na(transcript_length)) %>%
  group_by(screening_status) %>%
  summarise(
    count = n(),
    mean_length = round(mean(transcript_length)),
    median_length = round(median(transcript_length)),
    sd_length = round(sd(transcript_length)),
    min_length = min(transcript_length),
    max_length = max(transcript_length),
    .groups = 'drop'
  )

View(length_summary)
write.csv(length_summary, "non_screened_analyses/length_summary.csv", row.names = FALSE)

ggplot(complete_analysis %>% filter(!is.na(transcript_length)), 
       aes(x = screening_status, y = log10(transcript_length), fill = screening_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    title = "Transcript Length Comparison: Screened vs Non-Screened",
    x = "Screening Status",
    y = "Log10(Transcript Length)",
    fill = "Screening Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Visualization: Transcript-level MITF correlation comparison
ggplot(mitf_comparison %>% filter(!is.na(r_transcript)), 
                                       aes(x = screening_status, y = r_transcript, fill = screening_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    title = "Transcript-level MITF Correlation: Screened vs Non-Screened",
    x = "Screening Status",
    y = "Transcript Pearson Correlation with MITF",
    fill = "Screening Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Test GC content difference
if (nrow(gc_analysis) > 0 && sum(!is.na(complete_analysis$gene_gc_content)) > 0) {
  gc_screened <- complete_analysis %>% filter(screening_status == "CRISPR Screened", !is.na(gene_gc_content)) %>% pull(gene_gc_content)
  gc_nonscreened <- complete_analysis %>% filter(screening_status == "Non-CRISPR Screened", !is.na(gene_gc_content)) %>% pull(gene_gc_content)
  
  if (length(gc_screened) > 0 && length(gc_nonscreened) > 0) {
    gc_test <- wilcox.test(gc_screened, gc_nonscreened)
    cat("GC Content Wilcoxon test p-value:", format(gc_test$p.value, scientific = TRUE), "\n")
  }
}

# Test transcript length difference
length_screened <- complete_analysis %>% filter(screening_status == "CRISPR Screened", !is.na(transcript_length)) %>% pull(transcript_length)
length_nonscreened <- complete_analysis %>% filter(screening_status == "Non-CRISPR Screened", !is.na(transcript_length)) %>% pull(transcript_length)

if (length(length_screened) > 0 && length(length_nonscreened) > 0) {
  length_test <- wilcox.test(length_screened, length_nonscreened)
  cat("Transcript Length Wilcoxon test p-value:", format(length_test$p.value, scientific = TRUE), "\n")
}


####################### STATISTICAL TESTS
# Test GC content difference
if (nrow(gc_analysis) > 0 && sum(!is.na(complete_analysis$gene_gc_content)) > 0) {
  gc_screened <- complete_analysis %>% filter(screening_status == "CRISPR Screened", !is.na(gene_gc_content)) %>% pull(gene_gc_content)
  gc_nonscreened <- complete_analysis %>% filter(screening_status == "Non-CRISPR Screened", !is.na(gene_gc_content)) %>% pull(gene_gc_content)
  
  if (length(gc_screened) > 0 && length(gc_nonscreened) > 0) {
    gc_ttest <- t.test(gc_screened, gc_nonscreened)
    cat("GC Content t-test p-value:", format(gc_ttest$p.value, scientific = TRUE), "\n")
    cat("GC Content mean difference:", round(mean(gc_screened) - mean(gc_nonscreened), 4), "\n")
  }
}

# Test transcript length difference
length_screened <- complete_analysis %>% filter(screening_status == "CRISPR Screened", !is.na(transcript_length)) %>% pull(transcript_length)
length_nonscreened <- complete_analysis %>% filter(screening_status == "Non-CRISPR Screened", !is.na(transcript_length)) %>% pull(transcript_length)

if (length(length_screened) > 0 && length(length_nonscreened) > 0) {
  length_ttest <- t.test(length_screened, length_nonscreened)
  cat("Transcript Length t-test p-value:", format(length_ttest$p.value, scientific = TRUE), "\n")
  cat("Transcript Length mean difference:", round(mean(length_screened) - mean(length_nonscreened)), "bp\n")
}
