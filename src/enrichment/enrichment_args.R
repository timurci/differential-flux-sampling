parser = argparse::ArgumentParser(
    prog = "KEGG Enrichment Analysis"
)

parser$add_argument('zscore_path')

parser$add_argument('-k', '--key-column', required = T,
                    help = "name of the gene ID column")
parser$add_argument('-t', '--key-type', default = 'kegg',
                    choices = c('kegg', 'ncbi-geneid', 'ncbi-proteinid', 'uniprot'))
parser$add_argument('--sep', default = ";",
                    help = "separator character for multiple genes")

parser$add_argument('--organism', required = T,
                    help = "KEGG ID of the organism")

parser$add_argument('--repeats', default = 100, type = "integer",
                    help = "number of times to repeat tests (due to many number of ties)")
parser$add_argument('--repeat-threshold', default = .95, type = "double",
                    help = "percent threshold of repeats to include a pathway in list")

parser$add_argument('--p-cutoff', default = .05, type = "double")
parser$add_argument('--adjustment-method', default = "fdr",
                    choices = c('holm', 'hochberg', 'hommel', 'bonferroni',
                                'BH', 'BY', 'fdr', 'none'))

parser$add_argument('-o', '--output-dir', required = T)

cmd.args = NULL

do.debug = T
test.args = c('~/workspace/differential-flux-sampling/pyk.b10.r5/rank/zscore.csv',
              '--organism', 'ppu', '-k', 'gene_id', '--output-dir', 'pyk.b10.r5/enrichment',
              '--repeats', '10')

if (do.debug) {
    cmd.args = parser$parse_args(test.args)   
} else {
    cmd.args = parser$parse_args()
}

rm(parser, do.debug, test.args)

message(
    paste(
        format(names(cmd.args), width = 10),
        ":",
        cmd.args,
        collapse = "\n"
    )
)
