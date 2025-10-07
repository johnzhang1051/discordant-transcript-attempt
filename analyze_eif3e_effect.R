####################### Setting up environment
library(tidyverse)
library(biomaRt)
library(GenomicRanges)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

####################### LOAD DATA
avana_guide_map <- readr::read_csv("depmap-data/AvanaGuideMap.csv")
avana_logfold_change <- readr::read_csv("depmap-data/AvanaLogfoldChange.csv")
colnames(avana_logfold_change)[1] <- "sgRNA"
cell_info <- read.csv("depmap-data/Model.csv")
screen_sequence_map <- readr::read_csv("depmap-data/ScreenSequenceMap.csv")

# Define melanoma subtypes
melanoma_subtypes <- c("Melanoma", "Cutaneous Melanoma")

####################### GET EIF3E GUIDE EFFECTS

# Get EIF3E guides
eif3e_guides <- avana_guide_map %>%
  filter(Gene == "EIF3E (3646)")

# Get guide effects with cell line metadata
eif3e_effects <- avana_logfold_change %>%
  filter(sgRNA %in% eif3e_guides$sgRNA) %>%
  pivot_longer(cols = -sgRNA, names_to = "SequenceID", values_to = "log_fold_change") %>%
  left_join(screen_sequence_map, by = "SequenceID") %>%
  left_join(cell_info, by = "ModelID") %>%
  mutate(is_melanoma = OncotreeSubtype %in% melanoma_subtypes)

####################### MAP GUIDES TO TRANSCRIPTS

# Get EIF3E transcripts from Ensembl
ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
eif3e_transcripts <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name", "ensembl_exon_id",
                 "exon_chrom_start", "exon_chrom_end", "strand", "chromosome_name", "transcript_biotype"),
  filters = "external_gene_name",
  values = "EIF3E",
  mart = ensembl_data
)

# Prepare data for GenomicRanges overlap
eif3e_guides_map <- eif3e_guides %>%
  filter(grepl("^chr[0-9XY]", GenomeAlignment)) %>%
  separate(GenomeAlignment, into = c("chromosome", "guide_coordinate", "guide_strand"), sep = "_") %>%
  mutate(chromosome = gsub("chr", "", chromosome),
         guide_coordinate = as.integer(guide_coordinate),
         strand_char = guide_strand)

eif3e_exons <- eif3e_transcripts %>%
  mutate(strand_char = ifelse(strand == 1, "+", ifelse(strand == -1, "-", "*"))) %>%
  select(-strand)

# Find overlaps
exon_gr <- makeGRangesFromDataFrame(eif3e_exons, seqnames.field = "chromosome_name",
                                    start.field = "exon_chrom_start", end.field = "exon_chrom_end",
                                    strand.field = "strand_char", keep.extra.columns = TRUE)
guide_gr <- makeGRangesFromDataFrame(eif3e_guides_map, seqnames.field = "chromosome",
                                     start.field = "guide_coordinate", end.field = "guide_coordinate",
                                     strand.field = "strand_char", keep.extra.columns = TRUE)
overlap_hits <- findOverlaps(guide_gr, exon_gr, ignore.strand = TRUE)

eif3e_overlaps <- data.frame(
  guide_data = eif3e_guides_map[queryHits(overlap_hits), ],
  exon_data = eif3e_exons[subjectHits(overlap_hits), ]
)

####################### CLASSIFY GUIDES: DISCORDANT vs OTHER

# Load discordant list
discordant_list <- read.csv("resubmission_data/discordant_RESUBMISSION.csv")
colnames(discordant_list)[1] <- "transcript_id"

# Identify discordant EIF3E transcripts
eif3e_discordant <- unique(eif3e_overlaps$exon_data.ensembl_transcript_id)[
  unique(eif3e_overlaps$exon_data.ensembl_transcript_id) %in% discordant_list$transcript_id
]

# Classify guides
guides_discordant <- eif3e_overlaps %>%
  filter(exon_data.ensembl_transcript_id %in% eif3e_discordant) %>%
  pull(guide_data.sgRNA) %>% unique()

guides_other <- eif3e_overlaps %>%
  filter(!exon_data.ensembl_transcript_id %in% eif3e_discordant) %>%
  pull(guide_data.sgRNA) %>% unique()

# Add classification to effects
eif3e_guide_effects <- eif3e_effects %>%
  mutate(
    hits_discordant = sgRNA %in% guides_discordant,
    hits_other = sgRNA %in% guides_other,
    transcript_category = case_when(
      hits_discordant ~ "Discordant Targeted",
      hits_other ~ "Non-Discordant Only",
      TRUE ~ "Unclassified"
    )
  ) %>%
  filter(transcript_category != "Unclassified")

####################### GUIDE-LEVEL COMPARISON

# Aggregate by guide
guide_level_comparison <- eif3e_guide_effects %>%
  group_by(sgRNA, transcript_category, is_melanoma) %>%
  summarise(mean_lfc = mean(log_fold_change, na.rm = TRUE), 
            n_cell_lines = n(), .groups = "drop")

# Summary by category and cancer type
cancer_comparison <- guide_level_comparison %>%
  group_by(transcript_category, is_melanoma) %>%
  summarise(n_guides = n_distinct(sgRNA),
            mean_lfc = mean(mean_lfc, na.rm = TRUE),
            median_lfc = median(mean_lfc, na.rm = TRUE), .groups = "drop") %>%
  mutate(cancer_type = ifelse(is_melanoma, "Melanoma", "Non-Melanoma"))

print(cancer_comparison)

# T-tests for melanoma
if (any(guide_level_comparison$transcript_category == "Discordant Only" & guide_level_comparison$is_melanoma) &&
    any(guide_level_comparison$transcript_category == "Other Only" & guide_level_comparison$is_melanoma)) {
  
  disc_mel <- guide_level_comparison %>% filter(transcript_category == "Discordant Only", is_melanoma) %>% pull(mean_lfc)
  other_mel <- guide_level_comparison %>% filter(transcript_category == "Other Only", is_melanoma) %>% pull(mean_lfc)
  
  t_test_mel <- t.test(disc_mel, other_mel)
  cat("\nMelanoma T-test (Discordant Only vs Other Only):\n")
  cat("  P-value:", t_test_mel$p.value, "\n")
  cat("  Mean difference:", mean(disc_mel, na.rm = TRUE) - mean(other_mel, na.rm = TRUE), "\n")
}

# Visualization
ggplot(cancer_comparison, 
       aes(x = transcript_category, y = mean_lfc, fill = cancer_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "darkred", alpha = 0.5) +
  labs(title = "EIF3E Guide Effects: By Transcript Category",
       x = "Transcript Category", y = "Mean Log2 Fold Change", fill = "Cancer Type") +
  theme_minimal() +
  scale_fill_manual(values = c("Melanoma" = "red", "Non-Melanoma" = "gray70"))

# Export
write.csv(guide_level_comparison, "eif3e_analysis/guide_level_comparison.csv", row.names = FALSE)
write.csv(cancer_comparison, "eif3e_analysis/category_comparison.csv", row.names = FALSE)
