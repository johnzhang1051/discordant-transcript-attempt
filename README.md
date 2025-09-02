# discordant-transcript-attempt

This project is my first attempt at using Depmap data to identify discordant transcripts for MITF expression.

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
- AvanaGuide
  * Maps guide RNA's (Crispr) to genes and their locations 
[Ensembl](https://useast.ensembl.org/info/data/biomart/biomart_r_package.html)
- Biomart R Package
  * Lets us search genetic attributes (eg `exon_start_position`, `transcript_id`) by filtering to specific values such as gene, transcript, or exon id's


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

## Next Steps
- Analyze whether the guide targeting a transcript is correlated to the guide also affects the cancer cell
  * For this, I need to go through Depmap data to find Guide Effect data



## MITF High Low Analysis;
* In this side-work, we are looking at MITF-M expression distribution across Depmap cell-lines
* From there, we'll choose a cutoff for "high" and "low" expressing cell-lines
* We classify SequenceID's as high or low MITF expressing cell-lines
* Run separate transcript analyses for these groups, see if correlations are different
* I've chosen the median as the "cut-off" point
