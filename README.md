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

Biomart R Package
* Lets us search genetic attributes (eg `exon_start_position`, `transcript_id`) by filtering to specific values such as gene, transcript, or exon id's
* [Documentation](https://useast.ensembl.org/info/data/biomart/biomart_r_package.html)

## Analysis Flow:
1. `identify_MITF_high_low.R`:
   * Starts with Depmap cell-line data
   * Join in expression data for MITF-M (`ENST00000394351.9`)
   * Finds median expression for all cell-lines, using that as cutoff
   * Classifies Depmap cell-lines as MITF `high` or `low`
   * Exports list to `mitf_high_low/`
2. `analyze_sox10_effect.R`
   * Using Depmap guide-effect data and cell-lines
   * Find average guide-effect for melanoma and non-melanoma cell-lines, for sgRNA's that target SOX10
   * Verify that there is noticeable differences (melanoma effect should be much more negative) 
3. `calculate_guide_effects.R`
   * Starting with the discordant and correlated lists from the resubmission (`resubmission_data/`)
   * Map CRISPR guides from DepMap to transcript exons using GenomicRanges, count guides per transcript
   * Pull in guide effects (log2 fold change) from Depmap for each cell-line
   * Find average guide-effect for melanoma and non-melanoma cell-lines
   * Calculate FDR for guide-effect difference in melanoma vs. non-melanoma
   * Use results to find significant differences (`e.g FDR <= 0.05`)
4. `created_annotated_table.R`
   * Create "master" list of transcripts from correlated + discordant
   * Annotate the data + calculations we have
   * This will probably be used as the list of transcripts we do pooled screens on
5. `analyze_non_screened.R`
   * Compare Avana screened vs. non-screened transcripts
   * Specifically looking into the following:
      * Whether the transcript coordinates overlap (aren't unique)
      * Transcript length
      * % of GC content
   * Updates "master" table with new findings


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
- Compared screened vs. non-screened for GC %, Transcript Overlaps, and Transcript Length
- Analyzed SOX10 guides, found that guide effect was much stronger (negative) for melanoma cell-lines
- Using the resubmission discordant + correlated lists, I analyzed guide effects on Melanoma vs. Non-Melanoma cell-lines


## Next Steps
- Filter cell-lines to only MITF-M high expressing cell-lines, rerun guide effect analysis
- Create final annotated table of transcripts we want to screen
