# Microalgae-senoir-project
This is the Source code for the senior project, which contains

Bioinfirmatics Pipeline
- Genomics.py
  For running the Genomics Pipeline, which contains these tools
  QUAST -> BUSCO -> AUGUSTUS -> Extract Protein sequence for next tool -> DIAMOND -> Eggnog-mapper
  to extract the genome of all species.
- Transcriptomics.py
  For running the Transcriptomics Pipeline, the details of each tool are in the file "Transcriptomics_requirment"
- calculatetpm_all.sh
  To calculate the TPM (Transcriptome per Million) to use in WGCNA analysis, Labeling and Model Training.
- WGCNA_analysis.R
  For running a tool named WGCNA, which will use the TPM(Transcriptome per Million) to analyze the relationship between the gene and the important gene to calculate the score for define the label (High, Medium, Low)
