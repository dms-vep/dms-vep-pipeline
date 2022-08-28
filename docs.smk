"""``snakemake`` file that builds the HTML documentation for the pipeline.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import glob
import os
import re

# Variables and processing before building docs -------------------------------

# Get outputs of rules with `nb` as output (assumed Jupyter notebook) for docs.
rule_wildcards = {
    # for each rule with wildcards in nb output, define wildcards
    "fit_polyclonal": {
        "antibody_selection_group": antibody_selections["selection_group"].unique(),
    },
    "fit_globalepistasis": {
        "func_selection": func_selections["selection_name"].unique(),
    },
    "avg_antibody_escape": {
        "antibody": antibody_selections["antibody"].unique(),
    },
}
subindex_titles = {
    # name for each rule outputting wildcard notebooks
    "fit_polyclonal": "Fit ``polyclonal`` models",
    "fit_globalepistasis": "Fit global epistasis models",
    "avg_antibody_escape": "Antibody-escape averaged across replicates",
}
nbs = []
nblinks = {}
nbs_for_index = []
nb_subindices = {}
for name, ruleproxy in rules.__dict__.items():
    rule_nb = ruleproxy.output.get("nb")
    if rule_nb:
        if ruleproxy.rule.has_wildcards():
            assert set(rule_wildcards[name]) == set(ruleproxy.rule.wildcard_names)
            if len(rule_wildcards[name]) != 1:
                raise ValueError("currently only handles one wildcard per rule")
            else:
                wcs = list(rule_wildcards[name].values())[0]
                if len(wcs) == 0:
                    continue
            subindex = os.path.join(config["docs_source_dir"], f"{name}.rst")
            nb_subindices[subindex] = []
            for wc, nb in zip(wcs, expand(rule_nb, **rule_wildcards[name])):
                nbs.append(nb)
                nblink = os.path.join(config["docs_source_dir"], f"{name}_{wc}.nblink")
                nblinks[nblink] = nb
                nb_subindices[subindex].append(nblink)
        else:
            nbs.append(rule_nb)
            nblinks[os.path.join(config["docs_source_dir"], f"{name}.nblink")] = rule_nb
        nbs_for_index.append(name)

# If file should be linked, just specify file as value. If directory should be linked,
# specify value (files_in_directory, directory)
data_files = {
    "parental gene sequence": config["gene_sequence_codon"],
    "parental protein sequence": config["gene_sequence_protein"],
    "sequential-to-reference site numbers": config["site_numbering_map"],
    "codon-variant table": config["codon_variants"],
    **(
        {"processed barcode sequencing runs": config["processed_barcode_runs"]}
        if len(barcode_runs)
        else {}
    ),
    **(
        {"variant counts": (variant_count_files, config["variant_counts_dir"])}
        if variant_count_files
        else {}
    ),
}
if len(func_selections):
    # if we have functional selections, add those data files
    data_files.update(
        {
            "functional selection experiments": config["functional_selections"],
            "mutation effects for each functional selection": (
                muteffects_files,
                config["globalepistasis_dir"],
            ),
            "mutation functional effects replicate average (observed phenotype)": config[
                "muteffects_observed"
            ],
            "mutation functional effects replicate average (latent phenotype)": config[
                "muteffects_latent"
            ],
        }
    )
if len(antibody_selections):
    # if we have antibody selections, add those data files
    data_files.update(
        {
            "antibody selection experiments": config["antibody_selections"],
            "antibody escape values": (antibody_escape_files, config["escape_dir"]),
        }
    )
# add any extra data files specified in top-level Snakefile
data_files.update(extra_data_files)

# Rules ---------------------------------------------------------------------


rule make_graphs:
    """Build ``snakemake`` rulegraph and filegraph."""
    input:
        glob.glob("Snakefile*"),
        glob.glob("*.smk"),
        glob.glob(os.path.join(config["pipeline_path"], "Snakefile*")),
        glob.glob(os.path.join(config["pipeline_path"], "*.smk")),
    output:
        rulegraph=os.path.join(config["docs_source_dir"], "rulegraph.svg"),
        filegraph=os.path.join(config["docs_source_dir"], "filegraph.svg"),
    log:
        os.path.join(config["logdir"], "make_graphs.txt"),
    conda:
        "environment.yml"
    shell:
        """
        snakemake -F --rulegraph | dot -Tsvg > {output.rulegraph} 2> {log}
        snakemake -F --filegraph | dot -Tsvg > {output.filegraph} 2>> {log}
        """


rule make_nblink:
    """Make sphinx ``*.nblink`` file."""
    input:
        nb=lambda wc: nblinks[
            os.path.join(config["docs_source_dir"], f"{wc.nb}.nblink")
        ],
    output:
        nblink=os.path.join(config["docs_source_dir"], "{nb}.nblink"),
    params:
        nb_relpath=lambda _, input: os.path.relpath(
            input.nb, start=config["docs_source_dir"]
        ),
    log:
        os.path.join(config["logdir"], "make_nblink_{nb}.txt"),
    conda:
        "environment.yml"
    shell:
        """
        printf '{{\n    "path": "{params.nb_relpath}"\n}}' > {output.nblink} 2> {log}
        """


rule subindex:
    """Make ``*.rst`` subindex."""
    # regex "a^" never matches: https://stackoverflow.com/a/940840
    wildcard_constraints:
        subindex=(
            "|".join(re.escape(subindex) for subindex in nb_subindices)
        if nb_subindices
        else "a^"
        ),
    input:
        lambda wc: nb_subindices[wc.subindex],
    output:
        subindex="{subindex}",
    params:
        title=lambda _, output: subindex_titles[
            os.path.splitext(os.path.basename(output.subindex))[0]
        ],
    log:
        os.path.join(config["logdir"], os.path.basename("{subindex}") + ".txt"),
    conda:
        "environment.yml"
    script:
        "scripts/subindex.py"


rule docs_index:
    """Make ``index.rst`` file for sphinx docs."""
    input:
        [(f[0] if isinstance(f, tuple) else f) for f in data_files.values()],
        nbs,
        nb_subindices,
        nblinks,
        **(
            {
                "muteffects_observed": config["muteffects_observed_heatmap"],
                "muteffects_latent": config["muteffects_latent_heatmap"],
            }
            if len(func_selections)
            else {}
        ),
        rulegraph=rules.make_graphs.output.rulegraph,
        filegraph=rules.make_graphs.output.filegraph,
        antibody_escape_plots=antibody_escape_plots,
    output:
        index=os.path.join(config["docs_source_dir"], "index.rst"),
    params:
        docs_source_relpath=os.path.relpath(".", start=config["docs_source_dir"]),
        results_relpath=os.path.relpath(".", start=f"{config['docs']}/.."),
        data_files=data_files,
        nbs_for_index=nbs_for_index,
        have_func_selections=len(func_selections) > 0,
        have_antibody_selections=len(antibody_selections) > 0,
    log:
        os.path.join(config["logdir"], "docs_index.txt"),
    conda:
        "environment.yml"
    script:
        "scripts/docs_index.py"


rule sphinx_build:
    """Build sphinx docs."""
    input:
        rules.docs_index.input,
        index=rules.docs_index.output.index,
        conf=os.path.join(config["pipeline_path"], "conf.py"),
    output:
        docs=directory(config["docs"]),
    params:
        docs_source=lambda _, input: os.path.dirname(input.index),
        conf_path=lambda _, input: os.path.dirname(input.conf),
        project=config["description"],
        author=config["authors"],
        copyright=f"{config['authors']} ({config['year']})",
        # include HTML altair plots: https://stackoverflow.com/q/48889270
        html_extra_path=",".join(
            os.path.relpath(f, config["pipeline_path"])
            for f in [
                *(
                    [
                        config["muteffects_observed_heatmap"],
                        config["muteffects_latent_heatmap"],
                    ]
                    if len(func_selections)
                    else []
                ),
                *antibody_escape_plots,
            ]
        ),
    log:
        os.path.join(config["logdir"], "sphinx_build.txt"),
    conda:
        "environment.yml"
    shell:
        # The `-W` option means notebooks with wildcards that are in subindices must
        # have the "orphan" tag: https://nbsphinx.readthedocs.io/en/0.8.8/orphan.html
        """
        sphinx-build \
            -W \
            -b html \
            -c {params.conf_path} \
            -D project="{params.project}" \
            -D author="{params.author}" \
            -D copyright="{params.copyright}" \
            -D html_extra_path="{params.html_extra_path}" \
            {params.docs_source} \
            {output.docs} \
            &> {log}
        touch {output.docs}/.nojekyll
        """
