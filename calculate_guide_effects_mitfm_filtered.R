####################### Setting up environment
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(conflicted)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

####################### CONFIGURABLE INPUT - CHOOSE YOUR TRANSCRIPT LIST
# Modify these two lines to switch between different transcript lists:
transcript_list_file <- "resubmission_data/discordant_RESUBMISSION.csv"  # Changed to discordant
transcript_list_name <- "discordant_MITF_HIGH"  # Changed to discordant

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

gene_name_lookup <- results$structure %>%
  distinct(ensembl_transcript_id, external_gene_name)

# export Biomart data
#write.csv(results$structure, paste0("guide_effect/", transcript_list_name, "_exon_locations.csv"), row.names = FALSE)

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
    strand_char = guide_strand
  ) %>%
  dplyr::select(-guide_strand)  # Remove original strand column

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

guide_gr <- makeGRangesFromDataFrame(
  avana_guide_map,
  seqnames.field = "chromosome",
  start.field = "guide_coordinate",
  end.field = "guide_coordinate",
  strand.field = "strand_char",
  keep.extra.columns = TRUE
)

# Find overlaps
overlap_hits <- findOverlaps(guide_gr, exon_gr, ignore.strand = TRUE) #ignoring strand because shouldn't matter for Crispr

# Extract overlapping data
overlaps <- data.frame(
  guide_data = avana_guide_map[queryHits(overlap_hits), ],
  exon_data = exon_locations[subjectHits(overlap_hits), ]
)

# export overlaps data because it takes a while to run
#write.csv(overlaps, paste0("guide_effect/", transcript_list_name, "_exon_guide_overlaps.csv"), row.names = FALSE)

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
# Bring in MITF-M expression classifications, filter to "high-expressing"
mitf_classification <- read.csv("mitf_high_low/mitf_binary_classification.csv")
high_mitf_models <- mitf_classification %>%
  filter(mitf_binary == "High") %>%
  pull(ModelID)

# Define melanoma types
melanoma_subtypes <- c("Melanoma", "Cutaneous Melanoma")

# Filter to only guides that target our transcripts
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
  mutate(
    is_melanoma = OncotreeSubtype %in% melanoma_subtypes,
    # Only filter melanoma lines by MITF status
    keep_cell_line = case_when(
      is_melanoma ~ ModelID %in% high_mitf_models,  # Melanoma: only keep if MITF high
      !is_melanoma ~ TRUE                           # Non-melanoma: keep all
    )
  ) %>%
  filter(keep_cell_line) %>%  # Apply the filter
  select(-keep_cell_line)     # Remove helper column

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

###################### STATISTICAL TESTING FOR MELANOMA SPECIFICITY

# Perform t-tests for each transcript
transcript_stats <- transcript_guide_effects %>%
  group_by(exon_data.ensembl_transcript_id) %>%
  summarise(
    gene_name = first(exon_data.external_gene_name),
    n_guides = n_distinct(guide_data.sgRNA),
    n_melanoma_obs = sum(is_melanoma == TRUE),
    n_non_melanoma_obs = sum(is_melanoma == FALSE),
    
    # Calculate means
    mean_melanoma = mean(log_fold_change[is_melanoma == TRUE], na.rm = TRUE),
    mean_non_melanoma = mean(log_fold_change[is_melanoma == FALSE], na.rm = TRUE),
    median_melanoma = median(log_fold_change[is_melanoma == TRUE], na.rm = TRUE),
    median_non_melanoma = median(log_fold_change[is_melanoma == FALSE], na.rm = TRUE),
    
    # T-test
    t_test = list(tryCatch(
      t.test(
        log_fold_change[is_melanoma == TRUE],
        log_fold_change[is_melanoma == FALSE]
      ),
      error = function(e) list(p.value = NA, statistic = NA)
    )),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Extract t-test results
    p_value = sapply(t_test, function(x) x$p.value), # 
    effect_size = mean_melanoma - mean_non_melanoma
  ) %>%
  select(-t_test)

# Apply FDR correction
transcript_stats <- transcript_stats %>%
  mutate(
    fdr = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(effect_size) %>%
  rename(
    transcript_id = exon_data.ensembl_transcript_id)

###################### IDENTIFY INTERESTING TRANSCRIPTS

# Statistical approach: FDR-corrected significance + biological thresholds
interesting_transcripts <- transcript_stats %>%
  filter(
    fdr <= 0.1 &                          # Statistically significant (I made it 0.1 because there were some right on the edge)
    effect_size < 0 &                   # More essential in melanoma than non-melanoma
    mean_melanoma < 0                   # Actually essential in melanoma
  ) %>%
  arrange(fdr, effect_size) %>%
  select(transcript_id, gene_name, n_guides,
         mean_melanoma, mean_non_melanoma,
         effect_size, p_value, fdr)

print(interesting_transcripts)

transcript_stats <- transcript_stats %>%
  mutate(is_interesting = transcript_id %in% interesting_transcripts$transcript_id)

# Then when creating transcripts_without_guides:
transcripts_without_guides <- transcript_list %>%
  filter(!transcript_id_clean %in% transcripts_with_guides) %>%
  select(transcript_id = transcript_id_clean) %>%
  left_join(gene_name_lookup, by = c("transcript_id" = "ensembl_transcript_id")) %>%
  mutate(
    gene_name = external_gene_name,  # Use biomaRt gene name
    n_guides_targeting = 0,
    n_exons_targeted = 0,
    guide_list = ""
  ) %>%
  select(-external_gene_name)

# Combine to get complete transcript guide summary
transcript_guide_summary <- bind_rows(
  transcript_guide_counts,
  transcripts_without_guides
)

# Export results
combined_results <- transcript_guide_summary %>%
  left_join(
    transcript_stats %>% select(-gene_name),  # Remove duplicate gene_name from stats
    by = "transcript_id"
  ) %>%
  select(
    transcript_id,
    gene_name,
    n_guides,
    n_melanoma_obs,
    n_non_melanoma_obs,
    mean_melanoma,
    mean_non_melanoma,
    median_melanoma,
    median_non_melanoma,
    p_value,
    effect_size,
    fdr,
    is_interesting,
    n_exons_targeted,
    guide_list
  ) %>%
  # Fill NA for is_interesting with FALSE
  mutate(is_interesting = ifelse(is.na(is_interesting), FALSE, is_interesting)) %>%
  arrange(effect_size)
  

# Export combined results
write.csv(combined_results, 
          paste0("guide_effect/mitf_high/", transcript_list_name, "_guide_effect_results.csv"), 
          row.names = FALSE)

###################### Graphs

# Scatter plot with interesting transcripts highlighted
ggplot(transcript_stats,
       aes(x = mean_non_melanoma, 
           y = mean_melanoma,
           size = n_guides)) +
  geom_point(alpha = 0.5, shape = 21, fill = "gray80", color = "gray60", stroke = 0.3) +  # Non-interesting: gray
  geom_point(data = filter(transcript_stats, is_interesting),
             alpha = 0.7, shape = 21, 
             aes(fill = fdr),
             color = "black", stroke = 0.3) +  # Interesting: colored by FDR with black border
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "gray30", alpha = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray30", alpha = 0.5) +
  ggrepel::geom_text_repel(
    data = filter(transcript_stats, is_interesting),
    aes(label = transcript_id),
    size = 3,
    max.overlaps = 15,
    color = "black"
  ) +
  labs(
    title = paste(str_to_title(transcript_list_name), ": High MITF Expressing Melanoma vs Non-Melanoma"),
    subtitle = "Colored points with black borders = interesting targets; gray = others",
    x = "Guide Effect in Non-Melanoma",
    y = "Guide Effect in Melanoma",
    fill = "FDR",
    size = "Number of Guides"
  ) +
  theme_minimal() +
  scale_fill_gradient(low = "red", high = "yellow") +
  theme(legend.position = "right")

