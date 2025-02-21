# setwd('~/workspace/differential-flux-sampling/')

cmd.args = NULL
source('src/enrichment/enrichment_args.R')

if (!exists('cmd.args') || is.null(cmd.args)) quit(save = "no")

library(clusterProfiler)
library(tidyr)

scores = read.csv(cmd.args$zscore_path)

# 1. Threshold is applied to filter low-rank reactions.
# 2. Genes are grouped into "total", "positive", and "negative" based on z-score.
# 3. First occurrence (higher ranking) genes are selected in case of duplicates.
#
# Note: A gene can take part in the activity of multiple reactions.
#       It is possible to observe a positive flux shift in reaction X
#       and a negative flux shift in reaction Y associated with the same gene.

extract_genes = function(scores, key, threshold, sep = ";") {
  # Expand gene key column
  
  genes = separate_longer_delim(scores, key, delim = sep)
  genes = genes[genes[key] != '', ]
  
  genes
}

convert_genes = function(genes, key, method = "total") {
  # Filter genes by method and remove duplicates
  g = genes
  
  if (method == "positive") {
    g = subset(g, z.score > 0)
  } else if (method == "negative") {
    g = subset(g, z.score < 0)
  }
  
  g = g[!duplicated(g[,key]), ] # Keep the first occurrence
  
  sort(
    setNames(g$z.score, g[,key]),
    decreasing = T
  )
}

genes = extract_genes(scores, cmd.args$key_column, cmd.args$quantile, cmd.args$sep)

g.total = convert_genes(genes, cmd.args$key_column)
g.positive = convert_genes(genes, cmd.args$key_column, method = "positive")
g.negative = convert_genes(genes, cmd.args$key_column, method = "negative")


en.total.count = c()
en.positive.count = c()
en.negative.count = c()

update_count = function(counter, terms) {
  for (i in 1:length(terms)) {
    
    if (terms[i] %in% names(counter)) {
      counter[terms[i]] = counter[terms[i]] + 1
    } else {
      counter[terms[i]] = 1
    }
  }
  counter
}

for (i in 1:cmd.args$repeats) {
  
  en.total = gseKEGG(
    g.total,
    organism = cmd.args$organism,
    keyType = cmd.args$key_type,
    pvalueCutoff = cmd.args$p_cutoff,
    pAdjustMethod = cmd.args$adjustment_method,
    use_internal_data = FALSE,
  )
  
  en.total.count = update_count(en.total.count, en.total@result$ID)
  
  en.positive = gseKEGG(
      g.positive,
      organism = cmd.args$organism,
      keyType = cmd.args$key_type,
      pvalueCutoff = cmd.args$p_cutoff,
      pAdjustMethod = cmd.args$adjustment_method,
      use_internal_data = FALSE,
      scoreType = "pos",
  )
  
  en.positive.count = update_count(en.positive.count, en.positive@result$ID)
  
  en.negative = gseKEGG(
      g.negative,
      organism = cmd.args$organism,
      keyType = cmd.args$key_type,
      pvalueCutoff = cmd.args$p_cutoff,
      pAdjustMethod = cmd.args$adjustment_method,
      use_internal_data = FALSE,
      scoreType = "neg",
  )
  
  en.negative.count = update_count(en.negative.count, en.negative@result$ID)
  
}

en.total.select = names(en.total.count)[
  en.total.count >= max(en.total.count) * cmd.args$repeat_threshold
]
result.total = subset(en.total@result,
                      en.total.select %in% en.total@result$ID)

en.positive.select = names(en.positive.count)[
  en.positive.count >= max(en.positive.count) * cmd.args$repeat_threshold
]
result.positive = subset(en.positive@result,
                         en.positive.select %in% en.positive@result$ID)

en.negative.select = names(en.negative.count)[
  en.negative.count >= max(en.negative.count) * cmd.args$repeat_threshold
]
result.negative = subset(en.negative@result,
                         en.negative.select %in% en.negative@result$ID)

write.csv(result.total, file = sub('\\*', 'total', cmd.args$output))
write.csv(result.positive, file = sub('\\*', 'positive', cmd.args$output))
write.csv(result.negative, file = sub('\\*', 'negative', cmd.args$output))
