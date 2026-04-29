# Proteins Expressed in Predictive Ranges

This project contains R scripts designed to quantify and visualize how many proteins are expressed across different ranges of two predictive metrics: **TPM** (target RNA Transcripts Per Million) and **TRC** (Transcript Relative Consumption).

## Description

The analysis aims to understand the relationship between transcript-level predictors and actual protein expression. It classifies transcripts into bins based on their predicted values and then counts how many of those transcripts result in detectable protein expression.

Key features of the analysis:
- **Biomart Integration:** Uses `biomaRt` to filter for protein-coding transcripts and map Ensembl transcript IDs to UniProt accessions.
- **Quantile-based Binning:** Automatically divides the predictor values (TPC/TPM) into ranges based on quantiles.
- **Proteomics Integration:** Merges transcriptomic predictions with proteomics expression data.
- **Visualization:** Generates bar plots comparing the number of expressed proteins across the defined ranges for both predictors.

## Prerequisites

The scripts require R and the following libraries:
- `optparse`
- `biomaRt`
- `ggplot2`

You can install these packages using:
```R
install.packages(c("optparse", "ggplot2"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

## Project Structure

- `proteins_expressed_in_predictive_ranges.R`: Initial legacy version of the analysis script.
- `proteins_expressed_in_predictive_rangesv2.R`: Professional command-line version with parameterization, error handling, and modular functions.
- `plots/`: Directory containing generated visualization plots (PNG format).

## Usage

The `v2` script has been refactored into a professional command-line tool. You can now pass parameters directly from the terminal without editing the R script.

### Basic Execution

```bash
Rscript proteins_expressed_in_predictive_rangesv2.R \
  --trc_file path/to/trc_file.txt \
  --protein_file path/to/protein_file.txt
```

### Options

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-t`, `--trc_file` | Path to the TRC/TPM transcript file | **Required** |
| `-p`, `--protein_file` | Path to the Proteomics expression file | **Required** |
| `-o`, `--out_dir` | Output directory for plots | `plots/` |
| `--timepoint_prot` | Timepoint column name in the proteomics file | `UNTR_The_008_1` |
| `--timepoint_pred` | Timepoint identifier for the predictors | `UNTR_008_1` |
| `--normalize` | Normalize TRC and TPM values by their medians | `FALSE` |
| `--non_expressed` | Count non-expressed transcripts instead of expressed proteins | `FALSE` |
| `--disable_prot_dupl_filt`| Disable filtering of duplicate proteins | Filter enabled |
| `--disable_enst_dupl_filt`| Disable filtering of duplicate transcripts | Filter enabled |

### Example with advanced options

```bash
Rscript proteins_expressed_in_predictive_rangesv2.R \
  --trc_file data/UNTR_TRC_sample.txt \
  --protein_file data/Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt \
  --timepoint_prot UNTR_The_002_1 \
  --timepoint_pred UNTR_002_1 \
  --normalize \
  --out_dir my_custom_plots
```