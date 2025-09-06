####################### Setting up environment

library(dplyr)
MITF_gene <- "MITF (4286)" # we know this from the "clean_data.R" results
MITF_transcript <- "ENST00000394351.9" # we know this from the "clean_data.R" results

####################### Getting data
melanoma_transcript_expression <- readr::read_csv(
  file=file.path("cleaned_data", "melanoma_transcript_expression.csv")
)

melanoma_gene_expression <- readr::read_csv(
  file=file.path("cleaned_data", "melanoma_gene_expression.csv")
)
####################### Correlation analysis between MITF and transcripts:
transcript_names <- colnames(melanoma_transcript_expression)[!colnames(melanoma_transcript_expression) %in% c("model_id", "MITF", "index", "profile_id", "is_default_entry", MITF_transcript, "...1", "...2")]
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
    pearson_corr = cor(melanoma_transcript_expression[[MITF_transcript]], 
                       melanoma_transcript_expression[[transcript]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(melanoma_transcript_expression[[MITF_transcript]], 
                        melanoma_transcript_expression[[transcript]], 
                        method = "spearman", use = "complete.obs")
  )
})


# export correlations because they take forever to run
colnames(MITF_transcript_correlations)[1] <- "transcript_id"
write.csv(MITF_transcript_correlations,"correlation_results/MITF_transcript_correlations.csv", row.names = FALSE)

# also export correlations that are >0.5 for pearson and spearman
MITF_transcript_correlations_filtered <- MITF_transcript_correlations %>%
  filter(pearson_corr >= 0.5 & spearman_corr >= 0.5)
write.csv(MITF_transcript_correlations,"correlation_results/MITF_transcript_correlations_filtered.csv", row.names = FALSE)


####################### Correlation analysis for gene-level expression
gene_names <- colnames(melanoma_gene_expression)[!colnames(melanoma_gene_expression) %in% c("model_id", "MITF", "index", "profile_id", "is_default_entry", "...1", MITF_gene)]
total_genes <- length(gene_names)

MITF_gene_correlations <- map_dfr(gene_names, function(gene) {
  # Get current position
  current_pos <- which(gene_names == gene)
  
  # Print progress every 100 transcripts
  if(current_pos %% 100 == 0) {
    cat("Correlations finished/total:", current_pos, "/", total_genes, 
        "(", round(current_pos/total_genes*100, 1), "%)\n")
  }
  
  tibble(
    term = gene,
    pearson_corr = cor(melanoma_gene_expression[[MITF_gene]], 
                       melanoma_gene_expression[[gene]], 
                       method = "pearson", use = "complete.obs"),
    spearman_corr = cor(melanoma_gene_expression[[MITF_gene]], 
                        melanoma_gene_expression[[gene]], 
                        method = "spearman", use = "complete.obs")
  )
})

colnames(MITF_gene_correlations)[1] <- "gene_id"
# export correlations because they take forever to run
write.csv(MITF_gene_correlations,"correlation_results/MITF_gene_correlations.csv", row.names = FALSE)


#######################  Pull in correlations data (if working on it)
MITF_transcript_correlations <- readr::read_csv(
  file=file.path("correlation_results", "MITF_transcript_correlations.csv")
)
colnames(MITF_transcript_correlations)[1] <- "transcript_id"

MITF_gene_correlations <- readr::read_csv(
  file=file.path("correlation_results", "MITF_gene_correlations.csv")
)
colnames(MITF_gene_correlations)[1] <- "gene_id"

#######################  Map transcripts to parent genes

# using ensembl Biomart R package
# documentation: https://useast.ensembl.org/info/data/biomart/biomart_r_package.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
# Function to automatically query biomaRt with attribute grouping

ensembl_data <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

filtered_transcripts <- MITF_transcript_correlations %>%
  filter(pearson_corr > 0.5 | spearman_corr > 0.5)

transcript_ids <- sub("\\..*", "", sub("^[.]*", "", filtered_transcripts$transcript_id))



# Get the transcript-to-gene mapping data FIRST
transcript_gene_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl_data
)

# THEN process the mapping for CoCor
transcript_gene_map <- transcript_gene_mapping %>%
  dplyr::select(transcript_ID = ensembl_transcript_id, 
                Gene = external_gene_name) %>%
  dplyr::filter(!is.na(Gene)) %>%

####################### CoCor Analysis (ADDED)
library(cocor)
library(purrr)
library(tibble)

# Prepare data matrices for CoCor
trans_expr_raw <- melanoma_transcript_expression %>% 
  column_to_rownames("model_id")
trans_expr_raw <- trans_expr_raw[, !colnames(trans_expr_raw) %in% c("index", "profile_id", "is_default_entry", "...1", "...2")]

gene_expr_raw <- melanoma_gene_expression %>% 
  column_to_rownames("model_id")

# Find common samples and align both matrices
common_samples <- intersect(rownames(trans_expr_raw), rownames(gene_expr_raw))
trans_expr_raw <- trans_expr_raw[common_samples, ]
gene_expr_raw <- gene_expr_raw[common_samples, ]

# Reference: MITF transcript expression
y <- trans_expr_raw[[MITF_transcript]]

# Filter for valid pairs using the correct matrices
valid_pairs <- transcript_gene_map %>%
  rowwise() %>%
  filter(transcript_in_cols(transcript_ID, colnames(trans_expr_raw)),
         gene_in_cols(Gene, colnames(gene_expr_raw))) %>%
  ungroup() %>%
  filter(!is.na(Gene), Gene != "")

# Updated CoCor test function with length checking
run_cocor_test <- function(trans_id, gene_id, y_vec) {
  # Find columns in raw expression matrices
  trans_col <- colnames(trans_expr_raw)[grepl(trans_id, colnames(trans_expr_raw), fixed = TRUE)][1]
  gene_col <- colnames(gene_expr_raw)[grepl(gene_id, colnames(gene_expr_raw), fixed = TRUE)][1]
  
  if (is.na(trans_col) || is.na(gene_col)) {
    return(data.frame(
      transcript_ID = trans_id, Gene = gene_id, p_value = NA, n = 0,
      r_transcript = NA, r_gene = NA, r_transcript_gene = NA
    ))
  }
  
  # Extract vectors and ensure they're the same length
  x1 <- trans_expr_raw[[trans_col]]  # transcript expression
  x2 <- gene_expr_raw[[gene_col]]    # gene expression
  
  # Debug: check lengths
  if(length(x1) != length(x2) || length(x1) != length(y_vec)) {
    cat("Length mismatch - x1:", length(x1), "x2:", length(x2), "y:", length(y_vec), "\n")
    return(data.frame(
      transcript_ID = trans_id, Gene = gene_id, p_value = NA, n = 0,
      r_transcript = NA, r_gene = NA, r_transcript_gene = NA
    ))
  }
  
  complete_idx <- complete.cases(x1, x2, y_vec)
  x1 <- x1[complete_idx]
  x2 <- x2[complete_idx]
  y  <- y_vec[complete_idx]
  
  if (length(y) < 10) {
    return(data.frame(
      transcript_ID = trans_id, Gene = gene_id, p_value = NA, n = length(y),
      r_transcript = NA, r_gene = NA, r_transcript_gene = NA
    ))
  }
  
  # Rest of your function...
  r.jk <- cor(x1, y, method = "pearson")
  r.jh <- cor(x2, y, method = "pearson")
  r.kh <- cor(x1, x2, method = "pearson")
  
  test <- tryCatch({
    cocor.dep.groups.overlap(r.jk = r.jk, r.jh = r.jh, r.kh = r.kh, n = length(y))
  }, error = function(e) NULL)
  
  p_val <- tryCatch({
    if (!is.null(test)) as.numeric(test@steiger1980$p.value) else NA
  }, error = function(e) NA)
  
  return(data.frame(
    transcript_ID = trans_id, Gene = gene_id, p_value = p_val, n = length(y),
    r_transcript = r.jk, r_gene = r.jh, r_transcript_gene = r.kh
  ))
}

# Run CoCor analysis
cocor_results <- purrr::pmap_dfr(
  valid_pairs %>% dplyr::select(transcript_ID, Gene),
  function(transcript_ID, Gene) run_cocor_test(transcript_ID, Gene, y)
)

# Add FDR correction
cocor_results <- cocor_results %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

# Identify discordant transcripts
discordant_transcripts <- cocor_results %>%
  filter(!is.na(FDR) & FDR < 0.05,
         !is.na(r_transcript), !is.na(r_gene),
         r_transcript >= 0.5, r_gene < 0.5)


colnames(discordant_transcripts)[1] <- "transcript_id"
# Save results
write.csv(cocor_results, "correlation_results/cocor_results.csv", row.names = FALSE)
write.csv(discordant_transcripts, "correlation_results/discordant_transcripts.csv", row.names = FALSE)
