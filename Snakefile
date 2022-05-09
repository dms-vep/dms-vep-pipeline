"""``snakemake`` file that runs the pipeline.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import glob
import os
import textwrap

import pandas as pd

import helper_funcs

# Global variables and processing before pipeline -----------------------------

# Data frames for PacBio runs, Illumina barcode runs, antibody selections, etc.
# Some of these are written to CSV files, but only if they have changed.

pacbio_runs = helper_funcs.pacbio_runs_from_config(config["pacbio_runs"])

barcode_runs = helper_funcs.barcode_runs_from_config(
    config["barcode_runs"],
    valid_libraries=set(pacbio_runs["library"]),
)
os.makedirs(os.path.dirname(config["processed_barcode_runs"]), exist_ok=True)
helper_funcs.to_csv_if_changed(
    barcode_runs,
    config["processed_barcode_runs"],
    index=False,
)

antibody_selections = helper_funcs.get_antibody_selections(barcode_runs)
os.makedirs(os.path.dirname(config["antibody_selections"]), exist_ok=True)
helper_funcs.to_csv_if_changed(
    antibody_selections,
    config["antibody_selections"],
    index=False,
)

# Get BLAKE2b checksums and timestamps of *.csv files in `results`. Used
# below to re-adjust timestamps of some output files that haven't changed.
# Useful because some notebooks write output files for multiple samples
# only some of which may be changed from prior runs.
csv_times_checksums = {
    os.path.abspath(csv_file): {
        "checksum": helper_funcs.blake2b_checksum(csv_file),
        "ns": (os.stat(csv_file).st_atime_ns, os.stat(csv_file).st_mtime_ns),
    }
    for csv_file in glob.iglob("results/**/*.csv", recursive=True)
}

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


rule count_barcodes:
    """Count barcodes for a specific library-sample."""
    input:
        fastq_R1=(
            lambda wc: (
                barcode_runs.set_index("library_sample").at[
                    wc.library_sample, "fastq_R1"
                ]
            )
        ),
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


checkpoint variant_counts:
    """Get and analyze counts of different variants in each sample."""
    input:
        [
            os.path.join(config[f"barcode_{ftype}_dir"], f"{library_sample}.csv")
            for library_sample in barcode_runs["library_sample"]
            for ftype in ["counts", "counts_invalid", "fates"]
        ],
        config["gene_sequence_codon"],
        config["codon_variants"],
        config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/variant_counts.ipynb"),
    output:
        directory(config["variant_counts_dir"]),
        nb="results/notebooks/variant_counts.ipynb",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "variant_counts.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


def variant_count_files(wildcards):
    """Get variant count output files, and adjust timestamps of any not modified.

    Main goal is to back-modify timestamps of files created by `variant_counts` that were
    were not modified relative to start of pipeline. Somewhat hacky use of checkpointing.
    """
    subdir = checkpoints.variant_counts.get(**wildcards).output[0]
    count_files = {os.path.abspath(f) for f in glob.glob(f"{subdir}/*.csv")}
    expected_files = {
        os.path.abspath(f"{config['variant_counts_dir']}/{f}.csv")
        for f in barcode_runs.query("exclude_after_counts == 'no'")["library_sample"]
    }
    assert count_files == expected_files
    for f in count_files:
        if f in csv_times_checksums:
            if helper_funcs.blake2b_checksum(f) == csv_times_checksums[f]["checksum"]:
                os.utime(f, ns=csv_times_checksums[f]["ns"])
    return sorted(count_files)
