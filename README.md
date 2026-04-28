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
- `biomaRt`
- `ggplot2`

## Project Structure

- `proteins_expressed_in_predictive_ranges.R`: Initial version of the analysis script.
- `proteins_expressed_in_predictive_rangesv2.R`: Updated version with improved parameterization for different timepoints and support for analyzing non-expressed transcripts.
- `plots/`: Directory containing generated visualization plots (PNG format).

## Usage

The scripts are configured with several parameters at the top:

- `prot_dupl_filt`: Boolean to filter out duplicate proteins.
- `enst_dupl_filt`: Boolean to filter out duplicate transcripts.
- `normalize`: Boolean to normalize TRC and TPM values by their medians.
- `non_expressed`: (v2 only) If `TRUE`, the script counts transcripts that do *not* result in expressed proteins.
- `timepoint_prot` & `timepoint_pred`: (v2 only) Specify the identifiers for the proteomics and prediction samples respectively.

To run the analysis, ensure the file paths for the input TRC and protein tables are correctly set in the script, and then execute it in an R environment.