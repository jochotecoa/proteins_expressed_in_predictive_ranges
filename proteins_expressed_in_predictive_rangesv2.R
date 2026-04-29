#!/usr/bin/env Rscript

#' Proteins Expressed in Predictive Ranges
#'
#' This script quantifies and visualizes how many proteins are expressed 
#' across different ranges of predictive metrics (TPM and TRC).
#'
#' Usage: 
#' Rscript proteins_expressed_in_predictive_rangesv2.R --trc_file <path> --protein_file <path> ...

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(ggplot2)
})

#### Functions ####

#' Format data for ggplot
#' @param data Numeric vector of counts
#' @param fill String representing the predictor type (e.g., 'TRC', 'TPM')
#' @return A data frame formatted for ggplot
format_gg <- function(data, fill) {
  df <- as.data.frame(data)
  colnames(df) <- 'values'
  df$fill <- fill
  df$x <- rownames(df)
  return(df)
}

#' Create range columns based on quantiles
#' @param blocks Numeric vector of quantiles
#' @param predictor String column name of the predictor
#' @param table Data frame containing the predictor
#' @return Data frame with added boolean range columns
create_range_cols <- function(blocks, predictor, table) {
  for (i in seq_along(blocks)) {
    block_low_num <- blocks[i]
    colname <- paste0(predictor, '_higher_than_', block_low_num)
    block_group_log <- table[, predictor] >= block_low_num
    
    if (i != length(blocks)) {
      block_high_num <- blocks[i + 1]
      block_group2_log <- table[, predictor] < block_high_num
      block_group3_log <- block_group_log + block_group2_log
      block_groupfinal_log <- block_group3_log == 2
    } else {
      block_groupfinal_log <- block_group_log
    }
    table[, colname] <- block_groupfinal_log
  }
  return(table)
}

#### Main execution block ####

main <- function() {
  option_list <- list(
    make_option(c("-t", "--trc_file"), type = "character", default = NULL,
                help = "Path to the TRC/TPM transcript file", metavar = "character"),
    make_option(c("-p", "--protein_file"), type = "character", default = NULL,
                help = "Path to the Proteomics expression file", metavar = "character"),
    make_option(c("-o", "--out_dir"), type = "character", default = "plots",
                help = "Output directory for plots [default= %default]", metavar = "character"),
    make_option(c("--timepoint_prot"), type = "character", default = "UNTR_The_008_1",
                help = "Timepoint column name in the proteomics file [default= %default]", metavar = "character"),
    make_option(c("--timepoint_pred"), type = "character", default = "UNTR_008_1",
                help = "Timepoint identifier for the predictors [default= %default]", metavar = "character"),
    make_option(c("--normalize"), action = "store_true", default = FALSE,
                help = "Normalize TRC and TPM values by their medians"),
    make_option(c("--non_expressed"), action = "store_true", default = FALSE,
                help = "Count non-expressed transcripts instead of expressed proteins"),
    make_option(c("--disable_prot_dupl_filt"), action = "store_false", dest = "prot_dupl_filt", default = TRUE,
                help = "Disable filtering of duplicate proteins"),
    make_option(c("--disable_enst_dupl_filt"), action = "store_false", dest = "enst_dupl_filt", default = TRUE,
                help = "Disable filtering of duplicate transcripts")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$trc_file) || is.null(opt$protein_file)) {
    print_help(opt_parser)
    stop("Both --trc_file and --protein_file must be provided.\n", call. = FALSE)
  }
  
  if (!file.exists(opt$trc_file)) stop(paste("File not found:", opt$trc_file))
  if (!file.exists(opt$protein_file)) stop(paste("File not found:", opt$protein_file))
  
  if (!dir.exists(opt$out_dir)) {
    dir.create(opt$out_dir, recursive = TRUE)
  }
  
  message("--> Reading TRC sample file...")
  trc_table <- read.table(file = opt$trc_file, stringsAsFactors = FALSE, header = TRUE)
  
  if(!("TRC" %in% colnames(trc_table)) || !("targetRNA_TPM" %in% colnames(trc_table))) {
    # Attempt to handle if it doesn't have a header but we expect TRC and targetRNA_TPM
    # In the original script header wasn't specified, but usually it is needed to access $TRC
    message("Warning: Make sure TRC and targetRNA_TPM columns exist.")
  }
  
  trc_table$ensembl_transcript_id <- rownames(trc_table)
  
  if (opt$normalize) {
    message("--> Normalizing data by medians...")
    norm_factor_TRC <- median(trc_table$TRC, na.rm = TRUE)
    norm_factor_TPM <- median(trc_table$targetRNA_TPM, na.rm = TRUE)
    trc_table$TRC <- trc_table$TRC / norm_factor_TRC
    trc_table$targetRNA_TPM <- trc_table$targetRNA_TPM / norm_factor_TPM
  }
  
  message("--> Fetching protein-coding transcript biotypes from Ensembl biomaRt...")
  mart_human <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  
  trc_table_transcript_biotype <- getBM(attributes = c('ensembl_transcript_id', 'transcript_biotype'),
                                        filters = 'ensembl_transcript_id',
                                        values = trc_table$ensembl_transcript_id, 
                                        mart = mart_human)
  
  trc_table_2 <- merge(trc_table, trc_table_transcript_biotype, by = 'ensembl_transcript_id')
  trc_table_2 <- trc_table_2[trc_table_2$transcript_biotype == 'protein_coding', ]
  
  if (length(unique(trc_table_2$transcript_biotype)) > 1) {
    stop('Error: the transcript biotype has not been successfully filtered to protein_coding only.')
  }
  
  message("--> Classifying rows based on TRC/TPM values (Quantiles)...")
  blocks <- quantile(trc_table_2$TRC, seq(0, 0.9, 0.05), na.rm = TRUE)
  
  trc_table_2 <- create_range_cols(blocks = blocks, predictor = 'TRC', table = trc_table_2)
  trc_table_2 <- create_range_cols(blocks = blocks, predictor = 'targetRNA_TPM', table = trc_table_2)
  
  message("--> Fetching UniProt IDs for transcripts...")
  enst_unip <- getBM(attributes = c('ensembl_transcript_id', 'uniprot_gn'),
                     filters = 'ensembl_transcript_id',
                     values = trc_table_2$ensembl_transcript_id, 
                     mart = mart_human)
  
  trc_table_3 <- merge(trc_table_2, enst_unip, by = 'ensembl_transcript_id')
  
  message("--> Processing Protein expression values...")
  protein_table <- read.table(opt$protein_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  # Format protein names correctly
  protein_table <- protein_table[!grepl(pattern = ':', protein_table$Row.Names), ]
  names_split <- strsplit(as.character(protein_table$Row.Names), '\\|')
  protein_table$uniprot_gn <- sapply(names_split, function(x) if(length(x) > 1) x[2] else x[1])
  
  if (!(opt$timepoint_prot %in% colnames(protein_table))) {
    stop(paste("Timepoint column", opt$timepoint_prot, "not found in proteomics file."))
  }
  
  protein_table_sub <- protein_table[, c(opt$timepoint_prot, 'uniprot_gn')]
  
  trc_table_4 <- merge(x = trc_table_3, y = protein_table_sub, by = 'uniprot_gn', all.x = TRUE)
  tp_expr_col <- paste0(opt$timepoint_prot, '_expressed')
  trc_table_4[, tp_expr_col] <- trc_table_4[, opt$timepoint_prot] > 0
  trc_table_4[, tp_expr_col][is.na(trc_table_4[, tp_expr_col])] <- FALSE
  
  message("--> Calculating protein expression across ranges...")
  calculate_prots_in_ranges <- function(predictor_prefix) {
    n_prots_vec <- c()
    for (col in grep(paste0(predictor_prefix, '_higher'), colnames(trc_table_4), value = TRUE)) {
      group <- trc_table_4[trc_table_4[, col], ]
      
      if (opt$non_expressed) {
        group <- group[!group[, tp_expr_col], ]
      } else {
        group <- group[group[, tp_expr_col], ]
      }
      
      if (opt$enst_dupl_filt) { 
        group <- group[!duplicated(group$ensembl_transcript_id), ]
      }
      
      if (opt$prot_dupl_filt) {
        group_ids <- unique(group$uniprot_gn)
      } else {
        group_ids <- group$uniprot_gn
      }
      
      n_prots_vec <- c(n_prots_vec, length(group_ids))
      names(n_prots_vec)[length(n_prots_vec)] <- gsub(paste0(predictor_prefix, '_'), '', col)
    }
    return(n_prots_vec)
  }
  
  tpm_n_prots <- calculate_prots_in_ranges('targetRNA_TPM')
  trc_n_prots <- calculate_prots_in_ranges('TRC')
  
  message("--> Generating plot...")
  trc_n_prots_df <- format_gg(data = trc_n_prots, fill = 'TRC')
  tpm_n_prots_df <- format_gg(data = tpm_n_prots, fill = 'targetRNA_TPM')
  
  n_prots <- rbind(trc_n_prots_df, tpm_n_prots_df)
  n_prots$x <- as.numeric(gsub('higher_than_', '', n_prots$x))
  
  output_filename <- file.path(opt$out_dir, paste0('n_expr_prots_', opt$timepoint_prot, '_in_ranges_exsfgpr_', opt$timepoint_pred, '.png'))
  
  if (opt$non_expressed) {
    colnames(n_prots) <- c('Number_of_untranslated_transcripts', 'Predictor_type', 'Ranges_of_expression')
    p <- ggplot(n_prots, aes(y = Number_of_untranslated_transcripts, x = Ranges_of_expression, fill = Predictor_type)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      labs(title = "Untranslated Transcripts in Predictive Ranges", x = "Expression Range Threshold", y = "Count")
  } else {
    colnames(n_prots) <- c('Number_of_expressed_proteins', 'Predictor_type', 'Ranges_of_expression')
    p <- ggplot(n_prots, aes(y = Number_of_expressed_proteins, x = Ranges_of_expression, fill = Predictor_type)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      labs(title = "Expressed Proteins in Predictive Ranges", x = "Expression Range Threshold", y = "Count")
  }
  
  ggsave(filename = output_filename, plot = p, width = 8, height = 6, dpi = 300)
  message(paste("--> Plot successfully saved to:", output_filename))
}

# Execute main function
main()