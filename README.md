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

## Next Steps
- Figure out how to map Guides to Transcripts by using the exon locations from the previous step, likely using `AvanaGuide.csv` to map to CRISPR guides and positios to target
- Analyze whether the guide targeting a transcript is correlated to the guide also affects the cancer cell
