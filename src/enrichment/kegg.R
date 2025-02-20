# setwd('~/workspace/putida_sampling/')

cmd.args = NULL
source('src/enrichment/enrichment_args.R')

if (!exists('cmd.args') || is.null(cmd.args)) quit(save = "no")

library(clusterProfiler)
library(tidyr)

scores = read.table(cmd.args$zscore_path)

extract_genes = function(scores, key, threshold, sep = ";") {
  q = quantile(scores$rank, threshold)
  high.scores = subset(scores, rank <= q)
  
  genes = separate_longer_delim(high.scores, key, delim = sep)
  genes = genes[genes[key] != '', ]
  
  genes
}

genes = extract_genes(scores, cmd.args$key_column, cmd.args$quantile, cmd.args$sep)

en = enrichKEGG(
    setNames(genes$rank, genes[, cmd.args$key_column]),
    organism = cmd.args$organism,
    keyType = cmd.args$key_type,
    pvalueCutoff = cmd.args$p_cutoff,
    pAdjustMethod = cmd.args$adjustment_method,
    qvalueCutoff = cmd.args$q_cutoff,
    use_internal_data = FALSE
)

suffix = paste0("_q", cmd.args$quantile)

write(en@result$Description,
      file = paste("signf_descriptions",
                   suffix, ".txt", sep = ""))
