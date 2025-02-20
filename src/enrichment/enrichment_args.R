parser = argparse::ArgumentParser(
    prog = "KEGG Enrichment Analysis"
)

parser$add_argument('zscore_path')

parser$add_argument('-k', '--key-column', required = T,
                    help = "Name of the gene ID column")
parser$add_argument('-t', '--key-type', default = 'kegg',
                    choices = c('kegg', 'ncbi-geneid', 'ncbi-proteinid', 'uniprot'))
parser$add_argument('--sep', default = ";",
                    help = "Separator character for multiple genes")

parser$add_argument('--organism', required = T,
                    help = "Kegg ID of the organism")

parser$add_argument('--quantile', default = .1, type = "double",
                    help = "Reaction rank threshold")

parser$add_argument('--p-cutoff', default = .05, type = "double")
parser$add_argument('--adjustment-method', default = "fdr",
                    choices = c('holm', 'hochberg', 'hommel', 'bonferroni',
                                'BH', 'BY', 'fdr', 'none'))
parser$add_argument('--q-cutoff', default = .2, type = "double",
                    help = "q-value cutoff on enrichment tests")

parser$add_argument('-o', '--output', required = T)

cmd.args = NULL

do.debug = T
test.args = c('qweqew', '--organism', 'ppu', '-k', 'gene.set', '--output', 'test.out')

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