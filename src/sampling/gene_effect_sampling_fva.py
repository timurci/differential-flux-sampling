import cobra as cb
from cobra.flux_analysis.variability import flux_variability_analysis as fva
import json

def _optimum_fractional_bounds(model: cb.Model, rxn_ids: list[str], lb=0., ub=1.):
    var = fva(model, rxn_ids)

    var['minimum'] = var['maximum'] * lb
    var['maximum'] = var['maximum'] * ub

    var['bounds'] = list(zip(var['minimum'], var['maximum']))

    return var['bounds']

def active_rxn_bounds(model: cb.Model, genes: list[str], lb_active: float):
    """Get reaction bounds as a fraction of max. flux at the model's optimum"""
    rxns = set()

    # Find GPR associated reactions
    for g_id in genes:
        gene = model.genes.get_by_id(g_id)
        for rxn in gene.reactions:
            rxns.add(rxn)

    # Separate reactions by reversibility
    rev_rxns = {r for r in rxns if r.reversibility}
    irr_rxns = rxns.difference(rev_rxns)

    print(
        "Ignoring reversible reactions : ",
        json.dumps([r.id for r in rev_rxns], indent=2),
        sep=""
    )
    
    return _optimum_fractional_bounds(model, [r.id for r in irr_rxns], lb=lb_active)

def _tuple_format(tup, fmt='%.2f'):
    return '(' + ', '.join((fmt % f) for f in tup) + ')'

def _load_and_adjust_model(args):
    # Load model
    model = cb.io.read_sbml_model(args.model_path)
    
    print("Loaded model", model.id)
    print("Solver:", repr(model.solver))
    print("Model medium:", json.dumps(model.medium, indent=2))

    # Adjust GPR associated reactions of target genes
    if args.knockout:
        for g_id in args.target_genes:
            gene = model.genes.get_by_id(g_id)
            gene.knock_out()
            print("Knocked out", gene.id)
    else:
        bounds = active_rxn_bounds(model, args.target_genes, args.lb_target_rxns)
        for rxn_id, b in bounds.items():
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.bounds = b
            print("Set {:5} reaction bounds to".format(rxn.id),
                 _tuple_format(rxn.bounds))

    # Adjust biomass lower bound
    # Note: biomass is adjusted last to account for previous adjustments
    biomass = model.reactions.get_by_id(args.id_biomass)
    bounds = _optimum_fractional_bounds(model,
                                        args.id_biomass,
                                        lb=args.lb_biomass)
    biomass.bounds = bounds[biomass.id]
    print("Biomass reaction ID:", args.id_biomass)
    print("Set biomass reaction bounds to",
         _tuple_format(biomass.bounds))

    return model

def _parse_cmd_arguments():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        prog="gene effect by sampling",
        description="perform flux sampling on a GEM by activation or knockout of a set of genes",
    )

    parser.add_argument('model_path')
    parser.add_argument('-b', '--id-biomass',
                       help="biomass reaction by ID")
    parser.add_argument('--lb-biomass', default=.9, type=float,
                       help="lower bound of biomass by fraction of its optimum")
    parser.add_argument('--target-genes', nargs='+',
                       help="IDs of target genes to activate (or to knockout)")
    parser.add_argument('--lb-target-rxns', default=.05, type=float,
                       help="lower bound of target reactions by percentage of their optimum")
    parser.add_argument('--knockout', action='store_true',
                       help="knockout target genes")
    parser.add_argument('-n', '--samples', type=int,
                       help="Number of samples")
    parser.add_argument('-t', '--thinning', default=100, type=int,)
    parser.add_argument('--method', default="optgp", choices=["optgp", "achr"])
    parser.add_argument('-p', '--processes', default=1, type=int,
                       help="samples should be a multiple of number of processes")
    parser.add_argument('-o', '--output')

    args = parser.parse_args()
    model = _load_and_adjust_model(args)
    return args, model

def main():
    import logging
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    args, model = _parse_cmd_arguments()
    logger.debug("arguments: " + json.dumps(args.__dict__, indent=2))

    from time import time
    from datetime import timedelta
    from datetime import datetime

    start = time()
    logger.info("sampling starting at " + str(datetime.now()))

    samples = cb.sampling.sampling.sample(model,
                                          args.samples,
                                          args.method,
                                          args.thinning,
                                          args.processes)

    end = time()
    logger.info("sampling ended at " + str(datetime.now()))
    logger.info("elapsed time " + str(timedelta(seconds=(end - start))))

    samples.to_parquet(args.output)
    logger.info("sampling results are written into " + args.output)

if __name__ == "__main__":
    main()
