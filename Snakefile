configfile: "config.yaml"

def fformat(fraction):
    return "{:.0%}".format(fraction).replace('%', '')

PREFIX=".".join([
        config["scenario"],
        "b" + fformat(config["sampling"]["lb-biomass"]),
        "r" + fformat(config["sampling"]["lb-target"])
])

PATH_ENRICH="enrichment/" + ".".join([
        config["enrichment"]["organism"],
        "r" + str(config["enrichment"]["repeats"]),
        "t" + fformat(config["enrichment"]["repeat-threshold"])
])
 
rule all:
    input:
        f"{PREFIX}/sampling/WT.parquet",
        f"{PREFIX}/sampling/KO.parquet"

rule sampling:
    input:
        config["model"]["path"]
    output:
        f"{PREFIX}/" "sampling/{condition}.parquet"
    log:
        stdout=f"{PREFIX}/" "sampling/{condition}.stdout.log",
        stderr=f"{PREFIX}/" "sampling/{condition}.stderr.log"
    params:
        targets=expand("{gene}", gene=config["model"]["target-genes"])
    shell:
       " ".join([
            "python src/sampling/gene_effect_sampling_fva.py",
            "{input}",
            "-b {config[model][biomass]}",
            "--lb-biomass {config[sampling][lb-biomass]}",
            "--lb-target-rxns {config[sampling][lb-target]}",
            "--target-genes {params.targets}",
            "-n {config[sampling][samples]}",
            "-t {config[sampling][thinning]}",
            "-p {config[sampling][processes]}",
            "-o {output}",
            "> {log.stdout} 2> {log.stderr}"
        ])

rule rank:
    input:
        wt=f"{PREFIX}/" "sampling/WT.parquet",
        ko=f"{PREFIX}/" "sampling/KO.parquet"
    output:
        f"{PREFIX}/" "rank/distribution.png",
        csv=f"{PREFIX}/" "rank/zscore.csv"
    log:
        stdout=f"{PREFIX}/" "rank/zscore.stdout.log",
        stderr=f"{PREFIX}/" "rank/zscore.stderr.log"
    shell:
        " ".join([
            "python src/rank/zscore.py",
            "{config[model][path]}",
            "--cond1 {input.wt}",
            "--cond2 {input.ko}",
            "--gene-annotation {config[model][annotations]}",
            "-o {output.csv}",
            "> {log.stdout} 2> {log.stderr}"
        ])

rule enrichment:
    input:
        f"{PREFIX}/" "rank/zscore.csv"
    output:
        tot=f"{PREFIX}/" f"{PATH_ENRICH}.total.csv",
        pos=f"{PREFIX}/" f"{PATH_ENRICH}.positive.csv",
        neg=f"{PREFIX}/" f"{PATH_ENRICH}.negative.csv"
    log:
        stdout=f"{PREFIX}/" f"{PATH_ENRICH}.stdout.log",
        stderr=f"{PREFIX}/" f"{PATH_ENRICH}.stderr.log"
    shell:
        " ".join([
            "Rscript --vanilla src/enrichment/kegg.R",
            "{input}",
            "--key-column {config[enrichment][key]}",
            "--organism {config[enrichment][organism]}",
            "--repeats {config[enrichment][repeats]}",
            "--repeat-threshold {config[enrichment][repeat-threshold]}",
            "--output " f"{PREFIX}/" f"{PATH_ENRICH}.*.csv",
            "> {log.stdout} 2> {log.stderr}"
        ])
