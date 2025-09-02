####################### Setting up environment

library(depmap)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/Users/johnz/Documents/GitFiles/discordant-transcript-attempt")

####################### Getting data
# Get transcript/isoform expression data - so we know which transcripts are expressed
transcript_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv")
)
colnames(transcript_expression_data)[1] <- "index"

# Get cell line info - need this to identify melanoma cells (actually don't need this)
cell_info <- depmap_metadata()

# filter to all melanoma types that are not uveal, acral, mucosal, or uveal
specific_subtypes <- c("Melanoma", "Melanoma, amelanotic")
melanoma_cells <- cell_info %>%
  filter(subtype_disease %in% specific_subtypes)

####################### Filter data

# Get MITF expression data:
MITF_columns <- grep("ENST00000394351", colnames(transcript_expression_data), value = TRUE, ignore.case = TRUE)
MITF_expression_data <- transcript_expression_data %>%
  select(model_id, profile_id, all_of(MITF_columns))

####################### Combine data
combined_data <- melanoma_cells %>%
  inner_join(MITF_expression_data, by = c("depmap_id" = "model_id"))

combined_data <- combined_data %>%
  filter(!is.na(`ENST00000394351.9`))

####################### Generate distribution of MITF expression
# histogram for MITF expression
mean_val <- mean(combined_data$`ENST00000394351.9`, na.rm = TRUE)
median_val <- median(combined_data$`ENST00000394351.9`, na.rm = TRUE)
n_samples <- sum(!is.na(combined_data$`ENST00000394351.9`))

mitf_distribution_plot <- ggplot(combined_data, aes(x = `ENST00000394351.9`)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = mean_val, 
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = median_val, 
             color = "orange", linetype = "dashed", size = 1) +
  # Add text labels for mean and median
  annotate("text", x = mean_val, y = Inf, 
           label = paste("Mean =", round(mean_val, 2)), 
           color = "red", vjust = 2, hjust = -0.1, size = 3.5, fontface = "bold") +
  annotate("text", x = median_val, y = Inf, 
           label = paste("Median =", round(median_val, 2)), 
           color = "orange", vjust = 4, hjust = 1.1, size = 3.5, fontface = "bold") +
  # Add sample count label
  annotate("text", x = Inf, y = Inf, 
           label = paste("n =", n_samples), 
           color = "black", vjust = 2, hjust = 1.1, size = 3.5, fontface = "bold") +
  labs(
    title = "MITF Expression Distribution in Melanoma Cell Lines",
    subtitle = "Transcript: ENST00000394351.9",
    x = "MITF Expression (TPM Log+1)",
    y = "Number of Cell Lines",
    caption = "Red line: Mean, Orange line: Median"
  ) +
  theme_minimal()

print(mitf_distribution_plot)


# Based on results, I think we should use the median of 67.94 TPM Log + 1
# Standard deviation is quite high, wide range as well
library(psych)
describe(combined_data$ENST00000394351.9)



