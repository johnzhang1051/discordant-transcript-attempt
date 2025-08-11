####################### Installing depmap

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("depmap")


####################### Setting up environment

library(depmap)
library(tidyverse)
library(dplyr)
library(corrr)
library(ggplot2)
setwd("/Users/johnz/Documents/GitFiles/discordant-transcript-attempt")


####################### Getting data

# 1. Get cell line info - need this to identify melanoma cells:
cell_info <- depmap_metadata()
melanoma_cells <- filter(cell_info, grepl("Melanoma", subtype_disease, ignore.case = TRUE))

# 2. Get transcript/isoform expression data - so we know which transcripts are expressed
transcript_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv")
)
colnames(transcript_expression_data)[1] <- "index"

# 3. Get protein 
protein_gene_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
)
colnames(protein_gene_expression_data)[1] <- "model_id"


####################### Filter data

# Get MITF expression data:
MITF_columns <- grep("MITF", colnames(protein_gene_expression_data), value = TRUE, ignore.case = TRUE)
MITF_expression_data <- protein_gene_expression_data %>%
  select(model_id, all_of(MITF_columns))
colnames(MITF_expression_data)[2] <- "MITF"




####################### Combine data
combined_data <- MITF_expression_data %>%
  inner_join(transcript_expression_data, by = "model_id")

####################### Correlation analysis
# Calculate correlations between MITF and all transcripts
# Calculate correlations manually
transcript_names <- colnames(combined_data)[!colnames(combined_data) %in% c("model_id", "MITF", "index", "profile_id", "is_default_entry")]

total_transcripts <- length(transcript_names)

correlations <- map_dfr(transcript_names, function(transcript) {
  # Get current position
  current_pos <- which(transcript_names == transcript)
  
  # Print progress every 100 transcripts
  if(current_pos %% 100 == 0) {
    cat("Correlations finished/total:", current_pos, "/", total_transcripts, 
        "(", round(current_pos/total_transcripts*100, 1), "%)\n")
  }
  
  tibble(
    term = transcript,
    pearson_corr = cor(combined_data$MITF, combined_data[[transcript]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(combined_data$MITF, combined_data[[transcript]], 
                        method = "spearman", use = "complete.obs")
  )
})

head(spearman_corr)

# pearson_cors <- filtered_data %>%
#   select(-model_id, -index, -profile_id, -is_default_entry) %>%
#   correlate(method = "pearson") %>%
#   focus("MITF (4286)")
# 
# spearman_cors <- filtered_data %>%
#   select(-model_id, -index, -profile_id, -is_default_entry) %>%
#   correlate(method = "spearman") %>%
#   focus("MITF (4286)")

# Calculate gene-level correlations too
gene_correlations <- melanoma_genes %>%
  select(-ModelID) %>%
  correlate() %>%
  focus(MITF)

# Map transcripts to parent genes (simplified)
transcript_gene_map <- transcript_correlations %>%
  mutate(
    parent_gene = str_extract(term, "^[A-Z0-9]+"),  # Extract gene symbol
    transcript_corr = MITF
  ) %>%
  left_join(gene_correlations, by = c("parent_gene" = "term")) %>%
  rename(gene_corr = MITF.y) %>%
  select(-MITF.x)

# Find discordant transcripts (transcript correlates >20% better than gene)
discordant_transcripts <- transcript_gene_map %>%
  filter(
    abs(transcript_corr) > 0.5,  # Strong transcript correlation
    abs(transcript_corr) / abs(gene_corr) > 1.2  # 20% better than gene
  ) %>%
  arrange(desc(abs(transcript_corr)))

# See what you found
print(discordant_transcripts)


# Create scatter plot like their Figure 1A
ggplot(transcript_gene_map, aes(x = gene_corr, y = transcript_corr)) +
  geom_point(alpha = 0.6, color = "lightblue") +
  geom_point(data = discordant_transcripts, color = "red", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Transcript vs Gene Correlation with MITF",
    x = "Parent Gene Correlation with MITF", 
    y = "Transcript Correlation with MITF"
  ) +
  theme_minimal()
