####################### Guide Target Analysis and Cell Death Assessment
library(tidyverse)
library(ggplot2)
library(viridis)

# Set your transcript list name (change this as needed)
transcript_list_name <- "discordant_transcripts"
#transcript_list_name <- "correlated_transcripts"

####################### Load Results Data
# Load the guide effects results
guide_effects <- readr::read_csv(
  file = paste0("guide_effect/", transcript_list_name, "_guide_effects.csv")
)

# Load the detailed overlaps data
overlaps <- readr::read_csv(
  file = paste0("guide_effect/", transcript_list_name, "_exon_guide_overlaps.csv")
)

# Load original discordant transcript list
discordant_transcripts <- readr::read_csv(
  file=file.path("correlation_results", "discordant_transcripts.csv")
)
colnames(discordant_transcripts)[1] <- "transcript_id"

####################### 1. IDENTIFY TARGETED TRANSCRIPTS
# Basic targeting statistics
n_total_transcripts <- length(unique(guide_effects$exon_data.ensembl_transcript_id))
n_genes <- length(unique(guide_effects$gene_name))
total_guides <- sum(guide_effects$n_guides)

# Show transcripts with most/least guide coverage
cat("Transcripts with MOST guide coverage (top 10):\n")
top_coverage <- guide_effects %>% 
  arrange(desc(n_guides)) %>% 
  head(10) %>%
  select(gene_name, exon_data.ensembl_transcript_id, n_guides, mean_guide_effect)
print(top_coverage)


####################### 2. ASSESS CELL DEATH EFFECTS
# Define thresholds for effect classification
# I may change this into quartiles or something. Right now thresholds are kind of arbitray
essential_threshold <- -0.5    # Strong cell death
lethal_threshold <- -0.2       # Moderate cell death  
neutral_threshold <- 0.2       # Neutral
beneficial_threshold <- 0.5    # Potentially beneficial

# Classify transcripts by their effect
guide_effects <- guide_effects %>%
  mutate(
    effect_category = case_when(
      mean_guide_effect <= essential_threshold ~ "Highly Essential",
      mean_guide_effect <= lethal_threshold ~ "Moderately Essential", 
      mean_guide_effect <= neutral_threshold ~ "Neutral",
      mean_guide_effect <= beneficial_threshold ~ "Slightly Beneficial",
      TRUE ~ "Highly Beneficial"
    ),
    effect_category = factor(effect_category, levels = c(
      "Highly Essential", "Moderately Essential", "Neutral", 
      "Slightly Beneficial", "Highly Beneficial"
    ))
  )

# Summary statistics
effect_summary <- guide_effects %>%
  group_by(effect_category) %>%
  summarise(
    count = n(),
    percentage = round(100 * n() / nrow(guide_effects), 1),
    mean_effect = round(mean(mean_guide_effect), 3),
    .groups = "drop"
  )

cat("Effect Classification Summary:\n")
print(effect_summary)

# Key findings
essential_transcripts <- guide_effects %>% filter(mean_guide_effect <= lethal_threshold)
beneficial_transcripts <- guide_effects %>% filter(mean_guide_effect >= beneficial_threshold)

cat("- Essential transcripts (cause cell death):", nrow(essential_transcripts), 
    "(", round(100*nrow(essential_transcripts)/nrow(guide_effects), 1), "%)\n")
cat("- Beneficial transcripts (improve survival):", nrow(beneficial_transcripts),
    "(", round(100*nrow(beneficial_transcripts)/nrow(guide_effects), 1), "%)\n")

####################### 3. DETAILED TRANSCRIPT ANALYSIS

most_essential <- guide_effects %>%
  arrange(mean_guide_effect) %>%
  head(20) %>%
  select(gene_name, exon_data.ensembl_transcript_id, n_guides, mean_guide_effect, effect_category)

print(most_essential)


# Statistical test - are these transcripts more essential than average?
# One-sample t-test against neutral (0)
t_test_result <- t.test(guide_effects$mean_guide_effect, mu = 0)
cat("T-test p-value (vs neutral):", format(t_test_result$p.value, scientific = TRUE), "\n")
t_test_result$p.value
