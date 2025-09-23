# discordant-transcript-attempt

This project is my first attempt at using Depmap data to identify discordant transcripts for MITF expression, and finding interesting transcripts to study further.

For context, previously analyses have been done to find genes that correlated strongly to MITF expression.

However, due to alternative splicing and other processes, resulting RNA transcripts may correlate differently with MITF expression. This means if we ignore **genes** that don't correlate with MITF, we may **miss out on alternative transcripts from those genes that actually do correlate with MITF**.

Therefore, in this type of analysis we will try to find transcripts that are correlated with MITF, but come from genes that don't correlate with MITF.

This is based on previous work by Stephen Ostrowski.


## Data Sources:
[Depmap](https://depmap.org/portal/data_page/?tab=currentRelease)
- OmicsExpressionProteinCodingGenesTPMLogp1
  * How much of each gene is expressed in TPM units
- OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded
  * How much each transcript/isoform is expressed in TPM units
- OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded
  * How much each transcript/isoform is expressed in TPM units
- Model
   * Metadata describing all cancer models/cell lines within DepMap portal
- AvanaGuide
  * Maps guide RNA's (Crispr) to genes and their locations 
- ScreenSequenceMap
   * 
- Biomart R Package
  * Lets us search genetic attributes (eg `exon_start_position`, `transcript_id`) by filtering to specific values such as gene, transcript, or exon id's
  * [Documentation](https://useast.ensembl.org/info/data/biomart/biomart_r_package.html)

## Analysis Flow:
1. `clean_data.R`
   * Starts with Depmap transcript, gene expression, and cell-line data
   * Filter to only melanoma cell-lines for all 3 datasets
   * Clean column names
   * Identify MITF gene and MITF transcript column (important because version numbers change)
   * These are used in `identify_discordant.R`
2. `identify_MITF_high_low.R`:
   * Starts with Depmap cell-line data, filter to `Melanoma`, `Melanoma, amelanotic`
   * Join in expression data for MITF-M (`ENST00000394351.9`)
   * Finds median expression for all cell-lines, using that as cutoff
   * Classifies Depmap cell-lines as MITF `high` or `low`
   * Exports list to `mitf_high_low/`
3. `identify_discordant.R`
   * Starts with Depmap transcript and gene expression data
   * Using `CoCor` and `FDR < 0.05`, identifies discordant and correlated transcripts
   * Export to `correlation_results/`
   * NOTE: I'm thinking of using the MITF high/low and filtering to only MITF high for this code
      * But I'm not sure if we'll have enough data. I'll try it first to see.
4. `calculate_guide_effects.R`
   * Starting with discordant or correlated transcript list
   * Use BiomartR to get their exon coordinates
   * Import Depmap AvanaGuide Crispr data
   * Determine which transcripts were targeted in Crispr screen, and find `read_counts`
   * Export `guide_effect/...transcript_guide_counts.csv` to show # of guides targeting each transcript
   * From `read_counts`, calculate guide effect (see below for how to calculate)
   * Export guide_effect results
5. `analyze_guide_effects.R`
   * Starting with `guide_effect` data
   * List which transcripts had lowest guide_effect
   * We can conclude these are most "essential" for survival
   * Bc when these transcripts were screened, the cells died a lot
   * So likely these transcripts are important
   * And since we're looking at melanoma cell-lines, perhaps these transcripts are linked to melanoma
5. `analyze_MITF_high_low.R`
   * Starting with the `MITF_high_low` classifications
   * Reruns the whole Guide Effect analysis but keeps the MITF high/low classifications
   * Then reports whether there were differences in guide effects
6. `analyze_unique_promoters.R`
   * For list of transcripts, finds the # of promoters based on same data as resubmission
   * About ~60% of transcripts in our list don't have any promoter data on them
   * Exports `unique_promoters/all_transcripts_promoter_data.csv`, which is all the promoter data we have on a transcript-level
7. `analyze_non_screened.R`
   * Analyze transcripts that weren't covered in Avana screen
   * Specifically looking into the following:
      * Whether the transcript coordinates overlap (aren't unique)
      * PAM site availability
      * % of GC content
8. `created_annotated_table.R`
   * Create "master" list of transcripts from correlated + discordant
   * Annotate the data + calculations we have
   * This will probably be used as the list of transcripts we do pooled screens on


## Progress/Notes
- So far, I've gotten transcript-expression and protein-gene-expression data
- I've filtered the protein data to just MITF expression
- I've generated correlation values for MITF expression <-> transcript expression
- I've generated correlations between cell-line genes <-> MITF expression
- I've used Ensembl Biomart R package to map transcripts <-> parent genes <-> exons, and extracted exon locations for each transcript (note this part of the script takes a while)
- I created 2 R files:
  * `analysis_v1.R`: this is the old version that uses Ensembl's REST API to retrieve exon data (limited to 54k pulls, which is not enough)
  * `analysis_v2.R`: this is the newest version of my code that uses Ensembl's Biomart R package to pull all exon data instead (much more effective for the large amounts of data we need)
- I've used `GenomicRanges` library to find overlaps between `AvanaGuide` guide coordinates and the **transcript** coordinates, so we know which transcripts were targeted by the guides
- Drafted Guide effect calculation for transcripts targeted by guides (~81% discordant targeted by guides)
   * **Guide Effect is calculated by:**
   * Getting readcounts for each sample
   * Doing log2(read_count + 1) to get guide effect
   * For each cell-line, we find the median log2 value
   * Then we subtract that median from each value like so:
      * `(guide_effect_raw - guide_effect_median) = guide_effect_final`
   * Thus, the guide effect is normalized and centered 
- Determined expression cut-off for MITF-high vs. low expressing cell-lines
- Updated cell-line data to use updated from Depmap site
- Updated guide-effect to show # of guides targeting each transcript
- Identified # of promoters each transcript has

## Next Steps
* Recheck accuracy of data sources and calculations
* Analyzing non-screened transcripts (discordant and correlated)
   