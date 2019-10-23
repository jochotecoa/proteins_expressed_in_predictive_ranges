source('/share/script/hecatos/juantxo/analysis_trc/functions.R')
forceLibrary(c('biomaRt', 'ggplot2', 'dplyr'))


#### Functions ####
source(file = '/share/script/hecatos/juantxo/analysis_trc/functions.R')
format.gg = function(data, fill){
  # Formats a table ready to be ggplotted
  data = as.data.frame(data)
  colnames(data) = 'values'
  data$fill = fill
  data$x = rownames(data)
  return(data)
}
createRangeCols <- function(blocks, predictor, table) {
  for (i in 1:length(blocks)) {
    block_low.num = blocks[i]
    colname = paste0(predictor, '_higher_than_', block_low.num)
    block_group.log = table[, predictor] >= block_low.num
    if (i != length(blocks)) {
      block_high.num = blocks[i + 1]
      block_group2.log = table[, predictor] < block_high.num
      block_group3.log = block_group.log + block_group2.log
      block_groupfinal.log = block_group3.log == 2
    } else {
      block_groupfinal.log = block_group.log
    }
    table[, colname] = block_groupfinal.log
  }
  return(table)
}

#### Parameters ####

prot_dupl_filt = T
enst_dupl_filt = T
normalize = F
non_expressed = F
timepoint_prots = c('UNTR_The_002_1', 'UNTR_The_008_1', 'UNTR_The_024_1', 
                   'UNTR_The_072_1', 'UNTR_The_168_1', 'UNTR_The_240_1', 
                   'UNTR_The_336_1')
timepoint_pred = 'UNTR_002_1'

#### Pre-defined parameters ####

change_low = NULL
change_high = NULL

#### Get an example TRC sample file ####

setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/UNTR/TRCscore/')

files = list.files(pattern = timepoint_pred)
trc_table = read.table(file = files[1], stringsAsFactors = F) 
trc_table$ensembl_transcript_id = rownames(trc_table)
# Normalize data
if (normalize) {
  norm_factor_TRC = median(trc_table$TRC, na.rm = T)
  norm_factor_TPM = median(trc_table$targetRNA_TPM, na.rm = T)
  trc_table$TRC = trc_table$TRC / norm_factor_TRC
  trc_table$targetRNA_TPM = trc_table$targetRNA_TPM / norm_factor_TPM
}

#### Get only the protein-coding TRCs ####

mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'http://apr2018.archive.ensembl.org') 

# View(listAttributes(mart = mart.human))
trc_table_transcript_biotype = getBM(attributes = c('ensembl_transcript_id', 
                                                    'transcript_biotype'),
                           filters = 'ensembl_transcript_id',
                           values = trc_table$ensembl_transcript_id, 
                           mart = mart.human)

trc_table_2 = merge.data.frame(trc_table, trc_table_transcript_biotype, 
                               by = 'ensembl_transcript_id') %>%
  filter(transcript_biotype == 'protein_coding') %>%
  transcrToGene(aggregate = T)

# Check point
if (length(table(trc_table_2$transcript_biotype)) > 1) {
  message('Error: the transcript biotype has not been filtered')
}

#### Classify rows based on their TRC/TPM value ####

blocks = quantile(trc_table_2$targetRNA_TPM, seq(0,0.95,0.05))

trc_table_2 = trc_table_2 %>%
  createRangeCols(blocks = blocks, predictor = 'TRC', table = .) %>% 
  createRangeCols(blocks = blocks, predictor = 'targetRNA_TPM', table = .)


#### Get the proteins related to those rows ####

enst_unip = getBM(attributes = c('ensembl_gene_id',
                                 'uniprot_gn'),
                  filters = 'ensembl_gene_id',
                  values = trc_table_2$ensembl_gene_id, 
                  mart = mart.human)

trc_table_3 = merge.data.frame(trc_table_2, enst_unip, by = 'ensembl_gene_id')

#### Get the protein expression values ####
setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/')
setwd('Proteomics_Analyses_Cardiac_UNTR_GeneData/')

protein_table = read.table('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', 
                           header = T, sep = '\t')

# As always, format the protein names correctly
protein_table = protein_table %>% as_tibble() %>% 
  filter(!grepl(':', x = .$Row.Names))
protein_table$uniprot_gn = protein_table$Row.Names %>%
  as.character() %>% 
  strsplit('\\|') %>%
  lapply('[', 2) %>%
  as.character()

for (timepoint_prot in timepoint_prots) {
  # We only need our sample
  protein_table_0021 = protein_table[, c(timepoint_prot, 'uniprot_gn')]
  # Here, it is only all.x instead of all cause we don't want the proteins that 
  # are not in our table of study
  trc_table_4 = merge.data.frame(x = trc_table_3, y = protein_table_0021, 
                                 by = 'uniprot_gn', all.x = T)
  tp_expr_col = paste0(timepoint_prot, '_expressed')
  
  trc_table_4[, tp_expr_col] = trc_table_4[, timepoint_prot] > 0 
  trc_table_4[, tp_expr_col][is.na(trc_table_4[, tp_expr_col])] = F
  
  
  change_high_tmp_tpm = order(trc_table_4$targetRNA_TPM, decreasing = T) %>%
    .[1:100] %>%
    trc_table_4[., tp_expr_col] %>%
    sum()

  change_low_tmp_tpm = order(trc_table_4$targetRNA_TPM, decreasing = F) %>%
    .[1:100] %>%
    trc_table_4[, tp_expr_col][.] %>%
    sum()
  
  change_high_tmp_trc = order(trc_table_4$TRC, decreasing = T) %>%
    .[1:100] %>%
    trc_table_4[, tp_expr_col][.] %>%
    sum()
  
  change_low_tmp_trc = order(trc_table_4$TRC, decreasing = F) %>%
    .[1:100] %>%
    trc_table_4[, tp_expr_col][.] %>%
    sum()
  
  change_low_tmp = change_low_tmp_trc - change_low_tmp_tpm
  change_high_tmp = change_high_tmp_trc - change_high_tmp_tpm
  
  change_low = change_low %>% 
    rbind(c(change_low_tmp, timepoint_prot))
  colnames(change_low) = c('TRC_TPM_expressed_proteins_difference',
                           'timepoint_protein')
  
  change_high = change_high %>%
    rbind(c(change_high_tmp, timepoint_prot))
  colnames(change_high) = c('TRC_TPM_expressed_proteins_difference',
                           'timepoint_protein')
  
}

# rownames(change_high) = substr(x = rownames(change_high), 15, 17)
# rownames(change_low) = substr(x = rownames(change_low), 15, 17)

setwd("/share/script/hecatos/juantxo/analysis_trc/")
setwd('proteins_expressed_in_predictive_ranges/plots')

change_high = as.data.frame(change_high, stringsAsFactors = F)
colnames(change_high) = c('TRC_TPM_expressed_proteins_difference',
                          'timepoint_protein')
change_high[, 1] = as.numeric(change_high[, 1])

change_low = as.data.frame(change_low, stringsAsFactors = F)
colnames(change_low) = c('TRC_TPM_expressed_proteins_difference',
                          'timepoint_protein')
change_low[, 1] = as.numeric(change_low[, 1])

par(mfrow = c(1,2))
barplot(change_high[, 1], names.arg = change_high[, 2], 
        angle = 45, main = 'High end distribution differences', 
        xlab = 'Protein timepoint', 
        ylab = 'Difference in number of expressed proteins')
barplot(change_low$diff, names.arg = rownames(change_low), 
        angle = 45, main = 'Low end distribution differences', 
        xlab = 'Protein timepoint', 
        ylab = 'Difference in number of expressed proteins')
