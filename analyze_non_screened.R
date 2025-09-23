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
# Assuming you've already created this from your master table script
master_table <- read_csv("annotated_table/annotated_transcripts.csv")

# Filter for non-CRISPR screened transcripts
non_crispr_transcripts <- master_table %>%
  filter(crispr_screened == FALSE)

####################### Load existing biomaRt data from saved CSVs (no need to re-query!)
# Load your previously saved biomaRt exon location data
correlated_exon_locations <- read_csv("guide_effect/correlated_transcripts_exon_locations.csv")
discordant_exon_locations <- read_csv("guide_effect/discordant_transcripts_exon_locations.csv")

# Combine them
exon_locations <- bind_rows(
  correlated_exon_locations %>% mutate(source_list = "correlated"),
  discordant_exon_locations %>% mutate(source_list = "discordant")
)

# Clean transcript IDs (following your pattern)
non_crispr_transcripts$transcript_id_clean <- sub("\\..*", "", sub("^[.]*", "", non_crispr_transcripts$transcript_id))



####################### 1. COORDINATE OVERLAP ANALYSIS (exon-level overlaps)
cat("1. Analyzing coordinate overlaps at exon level...\n")

# Filter exon_locations to only non-CRISPR screened transcripts
exon_data <- exon_locations %>%
  filter(ensembl_transcript_id %in% non_crispr_transcripts$transcript_id_clean) %>%
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

# Get GC content data via biomaRt (since this wasn't in original queries)

# Set up biomaRt for GC content only
ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get unique transcript IDs for GC query
transcript_ids <- unique(non_crispr_transcripts$transcript_id_clean)

# Query for GC content only
gc_content_data <- tryCatch({
  getBM(
    attributes = c("ensembl_transcript_id", "percentage_gene_gc_content"),
    filters = "ensembl_transcript_id",
    values = transcript_ids,
    mart = ensembl_data
  )
}, error = function(e) {
  cat("Error getting GC content from biomaRt:", e$message, "\n")
  data.frame(ensembl_transcript_id = transcript_ids, percentage_gene_gc_content = NA)
})

# Use the GC content data we just retrieved
gc_analysis <- gc_content_data %>%
  rename(gene_gc_content = percentage_gene_gc_content) %>%
  filter(!is.na(gene_gc_content))

####################### Get GC content by extracting sequences directly
library(BSgenome.Hsapiens.UCSC.hg38)  # or hg19 if that's your reference
library(Biostrings)

# Filter to transcripts that have coordinate data
transcripts_with_coords <- transcript_coords_with_overlaps %>%
  filter(!is.na(chromosome_name), !is.na(transcript_start), !is.na(transcript_end))

if (nrow(transcripts_with_coords) > 0) {
  # Create GRanges for sequence extraction
  gr_for_gc <- GRanges(
    seqnames = paste0("chr", transcripts_with_coords$chromosome_name),  # Add 'chr' prefix for UCSC
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
      gene_gc_content = as.numeric(gc_content_values),  # Convert to percentage
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

####################### PAM SITE ANALYSIS (PLACEHOLDER FOR FUTURE)

# TODO: Implement PAM site analysis using DepMap guide data
# Potential approaches:
# - Use AvanaGuideMap.csv to analyze actual guide sequences and PAM sites
# - Calculate guide density per transcript from existing overlap data
# - Analyze guide efficiency scores if available in DepMap data

####################### COMBINE ALL ANALYSES
cat("4. Combining all technical analyses...\n")

# Merge all analyses
complete_analysis <- non_crispr_transcripts %>%
  # Join with coordinate/overlap analysis
  left_join(transcript_coords_with_overlaps, by = c("transcript_id_clean" = "ensembl_transcript_id")) %>%
  # Join with GC analysis
  left_join(gc_analysis, by = c("transcript_id_clean" = "ensembl_transcript_id"))


####################### VISUALIZATIONS

# Coordinate overlap distribution
complete_analysis %>%
  ggplot(aes(x = overlap_category, fill = transcript_classification)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Coordinate Overlap Patterns in Non-CRISPR Screened Transcripts",
    x = "Overlap Category", 
    y = "Number of Transcripts",
    fill = "Transcript Classification"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GC content distribution
if (nrow(gc_analysis) > 0 && sum(!is.na(complete_analysis$gene_gc_content)) > 0) {
  ggplot(complete_analysis %>% filter(!is.na(gene_gc_content)), 
               aes(x = gene_gc_content, fill = transcript_classification)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    geom_vline(xintercept = c(.20, .80), linetype = "dashed", color = "red") +
    labs(
      title = "GC Content Distribution in Non-CRISPR Screened Transcripts",
      x = "Gene GC Content (%)",
      y = "Number of Transcripts",
      fill = "Transcript Classification"
    ) +
    theme_minimal()
} else {
  cat("No GC content data available for plotting\n")
}

# Transcript length distribution by classification
complete_analysis %>%
  ggplot(aes(x = log10(transcript_length), fill = transcript_classification)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  labs(
    title = "Transcript Length Distribution in Non-CRISPR Screened Transcripts",
    x = "Log10(Transcript Length)",
    y = "Number of Transcripts", 
    fill = "Transcript Classification"
  ) +
  theme_minimal()
