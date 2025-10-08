####################### Setting up environment
library(tidyverse)
library(dplyr)
library(conflicted)

# Fix conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

####################### Load Required Data

# Load transcript lists
correlated_transcripts <- readr::read_csv(
  file = file.path("resubmission_data", "correlated_RESUBMISSION.csv")
)
colnames(correlated_transcripts)[1] <- "transcript_id"

discordant_transcripts <- readr::read_csv(
  file = file.path("resubmission_data", "discordant_RESUBMISSION.csv")
)
colnames(discordant_transcripts)[1] <- "transcript_id"

# Align column names
discordant_transcripts <- discordant_transcripts %>%
  rename(gene_name = Gene) %>%
  select(transcript_id, 
         gene_name, 
         p_value, 
         r_transcript, 
         r_gene, 
         FDR)

correlated_transcripts <- correlated_transcripts %>%
  rename(gene_name = Gene.x) %>%
  mutate(
    # r_transcript = max of all transcript correlations
    r_transcript = pmax(Pearson_Tsoi, Spearman_Tsoi, Pearson_CCLE, Spearman_CCLE, na.rm = TRUE),
    # r_gene = max of gene correlations
    r_gene = pmax(Pearson_gene, Spearman_gene, na.rm = TRUE)
    ) %>%
  select(transcript_id, 
         gene_name, 
         r_transcript, 
         r_gene)

# Load Guide effect data (MITF high)
discordant_guide_effects <- readr::read_csv(
  file = "guide_effect/mitf_high/discordant_MITF_HIGH_guide_effect_results.csv"
)
colnames(discordant_guide_effects)[1] <- "transcript_id"

correlated_guide_effects <- readr::read_csv(
  file = "guide_effect/mitf_high/correlated_MITF_HIGH_guide_effect_results.csv"
)
colnames(correlated_guide_effects)[1] <- "transcript_id"

####################### Create Master Table

# Add classification
correlated_transcripts <- correlated_transcripts %>%
  mutate(transcript_classification = "correlated")

discordant_transcripts <- discordant_transcripts %>%
  mutate(transcript_classification = "discordant")

# Combine all transcripts
all_transcripts <- bind_rows(correlated_transcripts, discordant_transcripts)

# Combine guide effects data
all_guide_effects <- bind_rows(
  discordant_guide_effects,
  correlated_guide_effects
) %>%
  rename(guide_effect_p_value = p_value) %>%
  rename(guide_effect_fdr = fdr)

# Create master table
master_table <- all_transcripts %>%
  left_join(all_guide_effects, by = "transcript_id") %>%
  mutate(crispr_screened = !is.na(n_guides)) %>%
  rename(gene_name = gene_name.x) %>%
  # Remove any remaining duplicate gene columns
  select(transcript_id, gene_name, transcript_classification, crispr_screened, r_transcript, r_gene, p_value, FDR, 
         n_guides, n_melanoma_obs, n_non_melanoma_obs, mean_melanoma, mean_non_melanoma, effect_size, guide_effect_p_value, guide_effect_fdr, is_interesting, n_exons_targeted, guide_list
  ) %>%
  arrange(desc(transcript_classification), desc(crispr_screened), effect_size, guide_effect_fdr)

####################### Save master table
write.csv(master_table, "annotated_table/annotated_transcripts.csv", row.names = FALSE)
