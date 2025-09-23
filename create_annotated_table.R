####################### Setting up environment
library(tidyverse)
library(dplyr)
library(conflicted)

# Fix select() conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

####################### Load Required Data
# Load transcript lists
correlated_transcripts <- readr::read_csv(
  file = file.path("correlation_results", "correlated_transcripts.csv")
)
colnames(correlated_transcripts)[1] <- "transcript_id"

discordant_transcripts <- readr::read_csv(
  file = file.path("correlation_results", "discordant_transcripts.csv")
)
colnames(discordant_transcripts)[1] <- "transcript_id"

# Load Guide targeting data
discordant_transcripts_guide_effects <- readr::read_csv(
  file = "guide_effect/discordant_transcripts_guide_effects.csv"
)
colnames(discordant_transcripts_guide_effects)[1] <- "transcript_id"

correlated_transcripts_guide_effects <- readr::read_csv(
  file = "guide_effect/correlated_transcripts_guide_effects.csv"
)
colnames(correlated_transcripts_guide_effects)[1] <- "transcript_id"

# Load Promoter data
all_transcripts_promoter_data <- readr::read_csv(
  file = "unique_promoters/all_transcripts_promoter_data.csv"
)
colnames(all_transcripts_promoter_data)[1] <- "transcript_id"

####################### Create Master Table

# Step 1: Create base transcript lists with classification
correlated_base <- correlated_transcripts %>%
  mutate(transcript_classification = "correlated")

discordant_base <- discordant_transcripts %>%
  mutate(transcript_classification = "discordant")

# Step 2: Combine all transcripts
all_transcripts <- bind_rows(correlated_base, discordant_base)

# Step 3: Combine guide effects data
all_guide_effects <- bind_rows(
  discordant_transcripts_guide_effects %>% mutate(source_classification = "discordant"),
  correlated_transcripts_guide_effects %>% mutate(source_classification = "correlated")
)

# Step 7: Create the master table
master_table <- all_transcripts %>%
  # Join with guide summary
  left_join(all_guide_effects, by = "transcript_id") %>%
  # Join with promoter summary
  left_join(all_transcripts_promoter_data, by = "transcript_id") %>%
  mutate(crispr_screened = !is.na(n_guides)) %>%
  # Clean up and organize columns
  select(
    transcript_id,
    transcript_classification,
    crispr_screened,
    n_promoters,
    everything(),
    -source_classification  # Remove if it exists from joins
  ) %>%
  # Sort by classification and transcript_id
  arrange(desc(transcript_classification), crispr_screened, desc(r_transcript))

####################### Save master table
write.csv(master_table, "annotated_table/annotated_transcripts.csv", row.names = FALSE)
