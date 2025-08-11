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


## Progress/Notes
- So far, I've gotten transcript-expression and protein-gene-expression data
- I've filtered the protein data to just MITF expression
- I've generated correlation values for MITF expression <-> transcript expression

## Next Steps
- Generate correlations between cell-line genes <-> MITF expression
- Identify discordant transcripts (where transcript<->MITF much stronger than gene<->MITF)
- Plot values
