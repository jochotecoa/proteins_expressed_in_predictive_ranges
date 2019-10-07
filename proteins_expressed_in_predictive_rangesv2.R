library('biomaRt')
library(ggplot2)

#### Functions ####
format.gg = function(data, fill){
  data = as.data.frame(data)
  colnames(data) = 'values'
  data$fill = fill
  data$x = rownames(data)
  return(data)
}
#### Parameters ####

prot_dupl_filt = T
enst_dupl_filt = T
normalize = F
non_expressed = T

#### Get an example TRC sample file ####

setwd('/share/analysis/hecatos/juantxo/Score/output/Output_Run_mrna_SEPT2019/')
setwd('V3/output/TRCscore/')

files = list.files(pattern = 'UNTR')
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
                               by = 'ensembl_transcript_id')
trc_table_2 = trc_table_2[trc_table_2$transcript_biotype == 'protein_coding', ]

# Check point
if (length(table(trc_table_2$transcript_biotype)) > 1) {
  message('Error: the transcript biotype has not been filtered')
}

#### Classify rows based on their TRC/TPM value ####

blocks = quantile(trc_table_2$TRC, seq(0,0.9,0.05))

createRangeCols <- function(blocks, predictor, table) {
  for (i in 1:length(blocks)) {
    block_low.num = blocks[i]
    colname = paste0(predictor, '_higher_than_', block_low.num)
    block_group.log = trc_table_2[, predictor] >= block_low.num
    if (i != length(blocks)) {
      block_high.num = blocks[i + 1]
      block_group2.log = trc_table_2[, predictor] < block_high.num
      block_group3.log = block_group.log + block_group2.log
      block_groupfinal.log = block_group3.log == 2
    } else {
      block_groupfinal.log = block_group.log
    }
    table[, colname] = block_groupfinal.log
  }
  return(table)
}

trc_table_2 = createRangeCols(blocks = blocks, 
                              predictor = 'TRC', 
                              table = trc_table_2)

trc_table_2 = createRangeCols(blocks = blocks, 
                              predictor = 'targetRNA_TPM', 
                              table = trc_table_2)

# for (i in 1:length(blocks)) {
#   block_low.num = blocks[i]
#   colname = paste0('TRC_higher_than_', block_low.num)
#   block_group.log = trc_table_2$TRC >= block_low.num
#   if (i != length(blocks)) {
#     block_high.num = blocks[i + 1]
#     block_group2.log = trc_table_2$TRC < block_high.num
#     block_group3.log = block_group.log + block_group2.log
#     block_groupfinal.log = block_group3.log == 2
#   } else {
#     block_groupfinal.log = block_group.log
#   }
#   trc_table_2[, colname] = block_groupfinal.log
# }
# 
# for (i in 1:length(blocks)) {
#   block_low.num = blocks[i]
#   colname = paste0('TPM_higher_than_', block_low.num)
#   block_group.log = trc_table_2$targetRNA_TPM >= block_low.num
#   if (i != length(blocks)) {
#     block_high.num = blocks[i + 1]
#     block_group2.log = trc_table_2$targetRNA_TPM < block_high.num
#     block_group3.log = block_group.log + block_group2.log
#     block_groupfinal.log = block_group3.log == 2
#   } else {
#     block_groupfinal.log = block_group.log
#   }
#   trc_table_2[, colname] = block_groupfinal.log
# }

#### Get the proteins related to those rows ####

enst_unip = getBM(attributes = c('ensembl_transcript_id',
                                 'uniprot_gn'),
      filters = 'ensembl_transcript_id',
      values = trc_table_2$ensembl_transcript_id, 
      mart = mart.human)

# dupl = duplicated(enst_unip$ensembl_transcript_id)
# dupl2 = enst_unip$ensembl_transcript_id[dupl]
# enst_unip2 = enst_unip[!(enst_unip$ensembl_transcript_id %in% dupl2), ]

trc_table_3 = merge.data.frame(trc_table_2, enst_unip, 
                               by = 'ensembl_transcript_id')

#### Get the protein expression values ####
setwd('/ngs-data/data/hecatos/Cardiac/Con_UNTR/Protein/')
setwd('Proteomics_Analyses_Cardiac_UNTR_GeneData/')

protein_table = read.table('Hecatos_Cardiac_Px_Untreated_pre-processed_renamed.txt', 
                           header = T, sep = '\t')

# As always, format the protein names correctly
protein_table = protein_table[!grepl(protein_table$Row.Names, pattern = ':'), ]
names = strsplit(as.character(protein_table$Row.Names), '\\|')
names = as.character(lapply(names, '[', 2))
protein_table$uniprot_gn = names
# We only need our sample
protein_table_0021 = protein_table[, c('UNTR_The_002_1', 'uniprot_gn')]
# Here, it is only all.x instead of all cause we don't want the proteins that 
# are not in our table of study
trc_table_4 = merge.data.frame(x = trc_table_3, y = protein_table_0021, 
                               by = 'uniprot_gn', all.x = T)
trc_table_4$UNTR_The_002_1_expr = trc_table_4$UNTR_The_002_1 > 0
trc_table_4$UNTR_The_002_1_expr[is.na(trc_table_4$UNTR_The_002_1_expr)] = F

#### Calculate how many proteins are expressed in those ranges ####
tpm_n_prots = NULL
i = 1
for (col in grep('TPM_higher', colnames(trc_table_4))) {
  
  tpmgroup = trc_table_4[trc_table_4[, col], ]
  if (non_expressed) {
    tpmgroup = tpmgroup[!tpmgroup$UNTR_The_002_1_expr, ]
  } else {
    tpmgroup = tpmgroup[tpmgroup$UNTR_The_002_1_expr, ]
  }
  # Filter out the transcripts that are associated to > 1 protein 
  if (enst_dupl_filt) { 
    tpmgroup = tpmgroup[!duplicated(tpmgroup$ensembl_transcript_id), ]
  }
  if (prot_dupl_filt) {
    tpmgroup_ids = unique(tpmgroup$uniprot_gn)
  } else {
    tpmgroup_ids = tpmgroup$uniprot_gn
  }
  tpmgroup_length = length(tpmgroup_ids)
  tpm_n_prots = c(tpm_n_prots, tpmgroup_length)
  names(tpm_n_prots)[length(tpm_n_prots)] = gsub('TPM_', '', 
                                                 colnames(trc_table_4)[col])
  i = i + 1
}

trc_n_prots = NULL
i = 1
for (col in grep('TRC_higher', colnames(trc_table_4))) {
  trcgroup = trc_table_4[trc_table_4[, col], ] # TRC in this range
  if (non_expressed) {
    trcgroup = trcgroup[!trcgroup$UNTR_The_002_1_expr, ]
  } else {
    trcgroup = trcgroup[trcgroup$UNTR_The_002_1_expr, ]
  }
  # Filter out the transcripts that are associated to > 1 protein 
  if (enst_dupl_filt) { 
    trcgroup = trcgroup[!duplicated(trcgroup$ensembl_transcript_id), ]
  }
  # A protein expressed by 2 transcripts is not 2 proteins
  if (prot_dupl_filt) {
    trcgroup_ids = unique(trcgroup$uniprot_gn)
  } else {
    trcgroup_ids = trcgroup$uniprot_gn
  }
  trcgroup_length = length(trcgroup_ids)
  trc_n_prots = c(trc_n_prots, trcgroup_length)
  names(trc_n_prots)[length(trc_n_prots)] = gsub('TRC_', '', 
                                                 colnames(trc_table_4)[col])
  i = i + 1
}

#### Plot data ####
trc_n_prots.df = format.gg(data = trc_n_prots, fill = 'TRC')
tpm_n_prots.df = format.gg(data = tpm_n_prots, fill = 'TPM')

n_prots = rbind(trc_n_prots.df, tpm_n_prots.df)
n_prots$x = as.numeric(gsub('higher_than_', '', n_prots$x))
if (non_expressed) {
  colnames(n_prots) = c('Number_of_untranslated_transcripts', 
                        'Predictor_type', 'Ranges_of_expression')
  ggplot(n_prots, aes(y = Number_of_untranslated_transcripts, 
                      x = Ranges_of_expression, fill = Predictor_type)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "Set1")
  
} else {
  colnames(n_prots) = c('Number_of_expressed_proteins', 
                        'Predictor_type', 'Ranges_of_expression')
  ggplot(n_prots, aes(y = Number_of_expressed_proteins, 
                      x = Ranges_of_expression, fill = Predictor_type)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "Set1")
  
}




