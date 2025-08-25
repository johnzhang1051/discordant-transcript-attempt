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

# 1. Get cell line info - need this to identify melanoma cells (actually don't need this)
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
protein_gene_expression_data <- protein_gene_expression_data %>%
  rename("MITF" = "MITF (4286)")

####################### Filter data

# Get MITF expression data:
MITF_columns <- grep("MITF", colnames(protein_gene_expression_data), value = TRUE, ignore.case = TRUE)
MITF_expression_data <- protein_gene_expression_data %>%
  select(model_id, all_of(MITF_columns))

####################### Combine data
combined_data <- MITF_expression_data %>%
  inner_join(transcript_expression_data, by = "model_id")

####################### Correlation analysis between MITF and transcripts
transcript_names <- colnames(combined_data)[!colnames(combined_data) %in% c("model_id", "MITF", "index", "profile_id", "is_default_entry")]

total_transcripts <- length(transcript_names)

MITF_transcript_correlations <- map_dfr(transcript_names, function(transcript) {
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

# export correlations because they take forever to run
write.csv(MITF_transcript_correlations,"MITF_transcript_correlations.csv", row.names = FALSE)

####################### Correlation analysis for gene-level expression
transcript_names <- colnames(protein_gene_expression_data)[!colnames(protein_gene_expression_data) %in% c("model_id")]

total_transcripts <- length(transcript_names)

MITF_gene_correlations <- map_dfr(transcript_names, function(transcript) {
  # Get current position
  current_pos <- which(transcript_names == transcript)
  
  # Print progress every 100 transcripts
  if(current_pos %% 100 == 0) {
    cat("Correlations finished/total:", current_pos, "/", total_transcripts, 
        "(", round(current_pos/total_transcripts*100, 1), "%)\n")
  }
  
  # tibble creates a small DF with the results (in this case the columns are: term, pearson_corr, spearman_corr)
  tibble(
    term = transcript,
    # Run correlations between MITF<->a transcript
    # complete.obs: removes pairs where any data is missing
    pearson_corr = cor(protein_gene_expression_data$MITF, protein_gene_expression_data[[transcript]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(protein_gene_expression_data$MITF, protein_gene_expression_data[[transcript]], 
                        method = "spearman", use = "complete.obs")
  )
})

# export correlations because they take forever to run
write.csv(MITF_gene_correlations,"MITF_gene_correlations.csv", row.names = FALSE)


#######################  Pull in correlations data (if working on it)
MITF_transcript_correlations <- readr::read_csv(
  file="MITF_transcript_correlations.csv"
)

colnames(MITF_transcript_correlations)[1] <- "transcript_id"

MITF_gene_correlations <- readr::read_csv(
  file="MITF_gene_correlations.csv"
)
colnames(MITF_gene_correlations)[1] <- "gene_id"

filtered_transcripts <- MITF_transcript_correlations %>%
  filter(pearson_corr > 0.5 | spearman_corr > 0.5)

#######################  Map transcripts to parent genes (not sure how to do this yet)

# using ensembl Rest API
# documentation: https://rest.ensembl.org/documentation/info/lookup

library(httr)
library(jsonlite)
library(xml2)
library(progress)

server <- "https://rest.ensembl.org"
all_ensembl_data <- data.frame()
failed_ids <- data.frame(transcript_id = character(), error = character(), stringsAsFactors = FALSE)

# Create progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed, ETA: :eta",
  total = nrow(filtered_transcripts),
  clear = FALSE,
  width = 60
)

loop_runs <- 0

for (x in filtered_transcripts$transcript_id) {
  loop_runs <- loop_runs + 1
  stable_transcript_id <- sub("\\..*", "", sub("^[.]*", "", x))
  ext <- paste("/lookup/id/", stable_transcript_id, "?expand=0", sep="")
  
  tryCatch({
    r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
    
    # Check specific status codes
    if (status_code(r) == 400) {
      failed_ids <- rbind(failed_ids, data.frame(
        transcript_id = stable_transcript_id, 
        error = "Bad Request (400)", 
        stringsAsFactors = FALSE
      ))
      pb$tick()
      next
    }
    
    if (status_code(r) == 404) {
      failed_ids <- rbind(failed_ids, data.frame(
        transcript_id = stable_transcript_id, 
        error = "Not Found (404)", 
        stringsAsFactors = FALSE
      ))
      pb$tick()
      next
    }
    
    # Will throw error if other HTTP error codes
    stop_for_status(r)
    
    ensembl_response_df <- data.frame(t(sapply(content(r), c)))
    
    # Add transcript_id for reference
    ensembl_response_df$transcript_id <- stable_transcript_id
    
    # Append to the main dataframe
    if (nrow(all_ensembl_data) == 0) {
      all_ensembl_data <- ensembl_response_df
    } else {
      # Use rbind.fill to handle column mismatches if needed
      all_ensembl_data <- rbind(all_ensembl_data, ensembl_response_df)
    }
    
  }, error = function(e) {
    # Catch any other errors (network issues, JSON parsing errors, etc.)
    failed_ids <<- rbind(failed_ids, data.frame(
      transcript_id = stable_transcript_id, 
      error = as.character(e$message), 
      stringsAsFactors = FALSE
    ))
  })
  
  pb$tick()
}

all_ensembl_data <- unique(all_ensembl_data)

#######################  Find discordant transcripts 
# (transcript correlates >20% better than gene)
# discordant_transcripts <- transcript_gene_map %>%
#   filter(
#     abs(transcript_corr) > 0.5,  # Strong transcript correlation
#     abs(transcript_corr) / abs(gene_corr) > 1.2  # 20% better than gene
#   ) %>%
#   arrange(desc(abs(transcript_corr)))
# 
# print(discordant_transcripts)

####################### Create scatter plot
# ggplot(transcript_gene_map, aes(x = gene_corr, y = transcript_corr)) +
#   geom_point(alpha = 0.6, color = "lightblue") +
#   geom_point(data = discordant_transcripts, color = "red", size = 2) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(
#     title = "Transcript vs Gene Correlation with MITF",
#     x = "Parent Gene Correlation with MITF", 
#     y = "Transcript Correlation with MITF"
#   ) +
#   theme_minimal()
