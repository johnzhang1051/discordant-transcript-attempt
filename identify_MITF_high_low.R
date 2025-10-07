####################### Setting up environment

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

####################### Getting data
# Get transcript/isoform expression data - so we know which transcripts are expressed
transcript_expression_data <- readr::read_csv(
  file=file.path("depmap-data", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv")
)
colnames(transcript_expression_data)[1] <- "index"

# Get cell line info - need this to identify melanoma cells (actually don't need this)
cell_info <- read.csv("depmap-data/Model.csv")


# filter to all melanoma types that are not uveal, acral, mucosal, or uveal
specific_subtypes <- c("Melanoma", "Melanoma, amelanotic", "Cutaneous Melanoma")
melanoma_cells <- cell_info %>%
  filter(OncotreeSubtype %in% specific_subtypes)

####################### Filter data

# Get MITF expression data:
MITF_columns <- grep("ENST00000394351", colnames(transcript_expression_data), value = TRUE, ignore.case = TRUE)
MITF_expression_data <- transcript_expression_data %>%
  select(model_id, profile_id, all_of(MITF_columns))

####################### Combine data
combined_data <- melanoma_cells %>%
  inner_join(MITF_expression_data, by = c("ModelID" = "model_id"))

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



############### ############### ############### Export classification of model-id's that are MITFM high or low

# Set cutoff for whether it's high or low
#cutoff_threshold <- median(combined_data$`ENST00000394351.9`, na.rm = TRUE)
cutoff_threshold <- 100
####################### Create Classifications

# Binary classification (High/Low using median)
mitf_classifications <- combined_data %>%
  mutate(
    # Method 1: Simple binary classification using median
    mitf_binary = case_when(
      `ENST00000394351.9` >= cutoff_threshold ~ "High",
      `ENST00000394351.9` < cutoff_threshold ~ "Low",
      TRUE ~ "Unknown"
    ),

    # Add the actual expression value for reference
    mitf_expression = `ENST00000394351.9`
  ) %>%
  # Select relevant columns for export
  select(
    ModelID,
    CellLineName,
    StrippedCellLineName,
    OncotreeSubtype,
    mitf_expression,
    mitf_binary
  )

####################### Summary Statistics

# Binary classification counts
binary_summary <- mitf_classifications %>%
  count(mitf_binary) %>%
  mutate(percentage = round(100 * n / sum(n), 1))
print(binary_summary)

####################### Export Classifications

# Export full classification table
full_output_file <- "mitf_high_low/mitf_expression_classifications_full.csv"
write.csv(mitf_classifications, full_output_file, row.names = FALSE)

# Export simplified binary classification (most commonly used)
binary_output <- mitf_classifications %>%
  select(ModelID, CellLineName, mitf_expression, mitf_binary)
binary_output_file <- "mitf_high_low/mitf_binary_classification.csv"
write.csv(binary_output, binary_output_file, row.names = FALSE)

####################### Create Visualization with Classifications

# Enhanced histogram showing classifications
ggplot(mitf_classifications, aes(x = mitf_expression)) +
  geom_histogram(aes(fill = mitf_binary), bins = 20, alpha = 0.7, color = "black") +
  geom_vline(xintercept = cutoff_threshold, 
             color = "red", linetype = "dashed", size = 1.2) +
  scale_fill_manual(values = c("High" = "darkred", "Low" = "steelblue")) +
  labs(
    title = "MITF-M Expression Classification in Melanoma Cell Lines",
    subtitle = paste("Binary classification using median threshold =", round(cutoff_threshold, 2)),
    x = "MITF Expression (TPM Log+1)",
    y = "Number of Cell Lines",
    fill = "MITF Classification",
    caption = "Red dashed line: Median threshold for High/Low classification"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
