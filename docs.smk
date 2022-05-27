"""``snakemake`` file that builds the HTML documentation for the pipeline.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import glob
import os

# Variables and processing before building docs -------------------------------

# Get outputs of rules with `nb` as output (assumed Jupyter notebook) for docs.
nb_outputs = {}
for name, ruleproxy in rules.__dict__.items():
    if ruleproxy.output.get("nb"):
        nb_outputs[name] = ruleproxy.output.get("nb")  # later expand with glob_wildcards

data_files = {
    "parental gene sequence": config["gene_sequence_codon"],
    "parental protein sequence": config["gene_sequence_protein"],
    "sequential-to-reference site numbers": config["site_numbering_map"],
    "codon-variant table": config["codon_variants"],
    "processed barcode sequencing runs": config["processed_barcode_runs"],
    "variant counts": config["variant_counts_dir"],
    "antibody selection experiments (grouped)": config["antibody_selection_groups"],
    "prob escapes for antibody selections": config["prob_escape_dir"],
    **extra_data_files,
}

# Rules ---------------------------------------------------------------------


rule make_graphs:
    """Build ``snakemake`` rulegraph, filegraph, and dag."""
    input:
        glob.glob("Snakefile*"),
        glob.glob("*.smk"),
        glob.glob(os.path.join(config["pipeline_path"], "Snakefile*")),
        glob.glob(os.path.join(config["pipeline_path"], "*.smk")),
    output:
        rulegraph=os.path.join(config["docs_source_dir"], "rulegraph.svg"),
        filegraph=os.path.join(config["docs_source_dir"], "filegraph.svg"),
        dag=os.path.join(config["docs_source_dir"], "dag.svg"),
    log:
        os.path.join(config["logdir"], "make_graphs.txt"),
    conda:
        "environment.yml"
    shell:
        """
        snakemake -F --rulegraph | dot -Tsvg > {output.rulegraph} 2> {log}
        snakemake -F --filegraph | dot -Tsvg > {output.filegraph} 2>> {log}
        snakemake -F --dag | dot -Tsvg > {output.dag} 2>> {log}
        """


rule docs_index:
    """Make ``index.rst`` file for sphinx docs."""
    input:
        rulegraph=rules.make_graphs.output.rulegraph,
        filegraph=rules.make_graphs.output.filegraph,
        dag=rules.make_graphs.output.dag,
    output:
        index=os.path.join(config["docs_source_dir"], "index.rst"),
    params:
        docs_source_relpath=os.path.relpath(".", start=config["docs_source_dir"]),
    log:
        os.path.join(config["logdir"], "docs_index.txt"),
    conda:
        "environment.yml"
    script:
        "scripts/docs_index.py"


rule sphinx_build:
    """Build sphinx docs."""
    input:
        rules.make_graphs.output,
        index=rules.docs_index.output.index,
        conf=os.path.join(config["pipeline_path"], "conf.py"),
    output:
        docs=directory(config["docs"]),
    params:
        docs_source=lambda _, input: os.path.dirname(input.index),
        conf_path=lambda _, input: os.path.dirname(input.conf),
        project=config["description"],
        author=config["authors"],
        copyright=f"{config['authors']} ({config['year']})"
    log:
        os.path.join(config["logdir"], "sphinx_build.txt"),
    conda: "environment.yml"
    shell:
        """
        sphinx-build \
            -b html \
            -W \
            -c {params.conf_path} \
            -D project="{params.project}" \
            -D author="{params.author}" \
            -D copyright="{params.copyright}" \
            {params.docs_source} \
            {output.docs} \
            &> {log}
        """
