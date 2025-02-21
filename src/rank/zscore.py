import pandas as pd
import cobra as cb
import seaborn as sns
import os
import logging


def _add_rank(zscore):
    zscore['rank'] = zscore['z-score'].abs().rank(ascending=False)

def _add_annotation(zscore, model, annotations):
    zscore['rxn_name'] = ''
    zscore['gene_id'] = ''

    for annot in annotations:
        zscore[annot] = ''

    for index in zscore.index:
        rxn = model.reactions.get_by_id(index)
        rxn_gene_id = []

        rxn_annotations = {annot: [] for annot in annotations}

        for g in rxn.genes:
            rxn_gene_id.append(g.id)
            for annot in rxn_annotations:
                rxn_annotations[annot].append(g.annotation.get(annot, ''))

        zscore.loc[index, 'rxn_name'] = rxn.name
        zscore.loc[index, 'gene_id'] = ";".join(rxn_gene_id)
        for k, v in rxn_annotations.items():
            zscore.loc[index, k] = ";".join(v)

def _shift_zscore(cond1_path, cond2_path):
    wt = pd.read_parquet(cond1_path)
    pt = pd.read_parquet(cond2_path)

    zscore = pt.sub(wt, axis=1).apply(
        lambda x : x.mean() / x.std() if not (x == 0).all() else 0,
        axis=0
    )

    zscore = zscore.loc[zscore != 0]

    return pd.DataFrame(
        zscore.sort_values(ascending=False, key=abs),
        columns=['z-score']
        )

def _parse_cmd_arguments():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Rank flux distribution shifts",
    )

    parser.add_argument('model_path')
    parser.add_argument('--cond1', required=True, help="default condition")
    parser.add_argument('--cond2', required=True, help="perturbed condition")
    parser.add_argument('--gene-annotation', default=[], nargs='*', help="additional gene annotation field(s)")
    parser.add_argument('-o', '--output', required=True)

    return parser.parse_args()

def main():
    import logging
    import json
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    args = _parse_cmd_arguments()
    logger.debug("arguments: " + json.dumps(args.__dict__, indent=2))

    # Calculate z-score
    zscore = _shift_zscore(args.cond1, args.cond2)

    # Add gene annotations and remove entries without any GPR rules
    model = cb.io.read_sbml_model(args.model_path)

    _add_annotation(zscore, model, args.gene_annotation)
    zscore = zscore.loc[zscore['gene_id'] != '']

    _add_rank(zscore)

    zscore.to_csv(args.output, index_label="rxn_id")

    # Plot z-score distribution
    base_dir = os.path.dirname(args.output)

    ax = sns.histplot(data=zscore, x="z-score", element="step")
    ax.get_figure().savefig(os.path.join(base_dir, "distribution.png"))


if __name__ == "__main__":
    main()
