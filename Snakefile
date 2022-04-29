"""``snakemake`` file that runs the pipeline.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import os
import textwrap

import pandas as pd

# Global variables extracted from config --------------------------------------
# Data frame of PacBio runs
pacbio_runs = (
    pd.read_csv(config["pacbio_runs"], dtype=str)
    .assign(pacbioRun=lambda x: x["library"] + "_" + x["run"].astype(str))
    .set_index("pacbioRun")
)
assert pacbio_runs.index.nunique() == len(pacbio_runs)


# Rules ---------------------------------------------------------------------


rule gene_sequence:
    """Get sequence of gene from PacBio amplicon."""
    input:
        gb=config["pacbio_amplicon"],
    output:
        codon=config["gene_sequence_codon"],
        prot=config["gene_sequence_protein"],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "gene_sequence.txt"),
    script:
        "scripts/gene_sequence.py"


rule align_parse_PacBio_ccs:
    """Align and parse PacBio CCS FASTQ file."""
    input:
        fastq=lambda wc: pacbio_runs.at[wc.pacbioRun, "fastq"],
        amplicon=config["pacbio_amplicon"],
        specs=config["pacbio_amplicon_specs"],
    output:
        outdir=directory(os.path.join(config["process_ccs_dir"], "{pacbioRun}")),
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "align_parse_PacBio_ccs_{pacbioRun}.txt"),
    script:
        "scripts/align_parse_PacBio_ccs.py"


rule analyze_pacbio_ccs:
    """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
    input:
        expand(rules.align_parse_PacBio_ccs.output.outdir, pacbioRun=pacbio_runs.index),
        config["pacbio_amplicon"],
        config["pacbio_amplicon_specs"],
        nb=os.path.join(config["pipeline_path"], "notebooks/analyze_pacbio_ccs.ipynb"),
    output:
        config["aligned_ccs_file"],
        nb="results/notebooks/analyze_pacbio_ccs.ipynb",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "analyze_pacbio_ccs.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


rule build_pacbio_consensus:
    """Build PacBio consensus sequences for barcodes."""
    input:
        config["aligned_ccs_file"],
        config["gene_sequence_codon"],
        nb=os.path.join(
            config["pipeline_path"], "notebooks/build_pacbio_consensus.ipynb"
        ),
    output:
        config["nt_variants"],
        nb="results/notebooks/build_pacbio_consensus.ipynb",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "build_pacbio_consensus.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


rule build_codon_variants:
    """Build codon-variant table."""
    input:
        config["nt_variants"],
        config["gene_sequence_codon"],
        config["site_numbering_map"],
        config["mutation_design_classification"],
        config["neut_standard_barcodes"],
        nb=os.path.join(config["pipeline_path"], "notebooks/build_codon_variants.ipynb"),
    output:
        config["codon_variants"],
        nb="results/notebooks/build_codon_variants.ipynb",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "build_codon_variants.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


checkpoint process_barcode_runs:
    """Check barcode runs and write CSV with sample names and other information."""
    input:
        csv=config["barcode_runs"],
    output:
        csv=config["processed_barcode_runs"],
    params:
        pacbio_libraries=pacbio_runs["library"].tolist(),
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "process_barcode_runs.txt"),
    script:
        "scripts/process_barcode_runs.py"


def barcode_count_fate_csvs(wildcards):
    """Get list of all barcode count CSVs."""
    fname = checkpoints.process_barcode_runs.get(**wildcards).output.csv
    library_samples = pd.read_csv(fname)["library_sample"]
    return [
        os.path.join(config[f"barcode_{ftype}_dir"], f"{f.csv}")
        for f in library_samples
        for ftype in ["counts", "counts_invalid", "fates"]
    ]


def barcode_fastq_R1(wildcards):
    """Get R1 FASTQs for a specific library-sample."""
    fname = checkpoints.process_barcode_runs.get(**wildcards).output.csv
    return [
        f.strip()
        for f in pd.read_csv(fname)
        .set_index("library_sample")
        .at[wildcards.library_sample, "fastq_R1"]
        .split(";")
    ]


rule count_barcodes:
    """Count barcodes for a specific library-sample."""
    input:
        fastq_R1=barcode_fastq_R1,
        variants=config["codon_variants"],
    output:
        counts=os.path.join(config["barcode_counts_dir"], "{library_sample}.csv"),
        counts_invalid=os.path.join(
            config["barcode_counts_invalid_dir"], "{library_sample}.csv"
        ),
        fates=os.path.join(config["barcode_fates_dir"], "{library_sample}.csv"),
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "count_barcodes_{library_sample}.txt"),
    script:
        "scripts/count_barcodes.py"
