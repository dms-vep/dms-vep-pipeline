"""``snakemake`` file that includes the pipeline code.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import glob
import os


include: "funcs.smk"  # import functions


# Global variables and processing before pipeline -----------------------------

# Data frames for PacBio runs, Illumina barcode runs, antibody selections, etc.
# Some of these are written to CSV files, but only if they have changed.

pacbio_runs = pacbio_runs_from_config(config["pacbio_runs"])

barcode_runs = barcode_runs_from_config(
    config["barcode_runs"],
    valid_libraries=set(pacbio_runs["library"]),
)
os.makedirs(os.path.dirname(config["processed_barcode_runs"]), exist_ok=True)
to_csv_if_changed(barcode_runs, config["processed_barcode_runs"], index=False)

library_sample_to_library = barcode_runs.set_index("library_sample")[
    "library"
].to_dict()
library_sample_to_sample = barcode_runs.set_index("library_sample")["sample"].to_dict()

variant_count_files = [
    os.path.join(config["variant_counts_dir"], f"{library_sample}.csv")
    for library_sample in barcode_runs.query("exclude_after_counts == 'no'")[
        "library_sample"
    ]
]

antibody_selections = get_antibody_selections(barcode_runs)
os.makedirs(os.path.dirname(config["antibody_selections"]), exist_ok=True)
to_csv_if_changed(antibody_selections, config["antibody_selections"], index=False)

antibody_selection_group_samples = {
    selection_group: sorted(
        set(
            antibody_selections.query("selection_group == @selection_group")[
                ["antibody_library_sample", "no-antibody_library_sample"]
            ].values.ravel()
        )
    )
    for selection_group in antibody_selections["selection_group"].unique()
}

prob_escape_files = [
    os.path.join(config["prob_escape_dir"], f"{selection_group}_{suffix}.csv")
    for selection_group in antibody_selections["selection_group"].unique()
    for suffix in ["prob_escape", "neut_standard_fracs", "neutralization"]
]


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
    params:
        parser_params=config["illumina_barcode_parser_params"],
        library=lambda wc: barcode_runs.set_index("library_sample").at[
            wc.library_sample, "library"
        ],
        sample=lambda wc: barcode_runs.set_index("library_sample").at[
            wc.library_sample, "sample"
        ],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "count_barcodes_{library_sample}.txt"),
    script:
        "scripts/count_barcodes.py"


rule variant_counts:
    """Get counts of variants for each sample."""
    input:
        barcode_counts=rules.count_barcodes.output.counts,
        codon_variants=config["codon_variants"],
        gene_sequence_codon=config["gene_sequence_codon"],
    output:
        counts=os.path.join(config["variant_counts_dir"], "{library_sample}.csv"),
    params:
        library=lambda wc: barcode_runs.set_index("library_sample").at[
            wc.library_sample, "library"
        ],
        sample=lambda wc: barcode_runs.set_index("library_sample").at[
            wc.library_sample, "sample"
        ],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "variant_counts_{library_sample}.txt"),
    script:
        "scripts/variant_counts.py"


rule analyze_variant_counts:
    """Analyze counts of different variants in each sample."""
    input:
        expand(
            rules.count_barcodes.output.counts,
            library_sample=barcode_runs["library_sample"],
        ),
        expand(
            rules.count_barcodes.output.counts_invalid,
            library_sample=barcode_runs["library_sample"],
        ),
        expand(
            rules.count_barcodes.output.fates,
            library_sample=barcode_runs["library_sample"],
        ),
        variant_count_files,
        config["gene_sequence_codon"],
        config["codon_variants"],
        config["site_numbering_map"],
        config["processed_barcode_runs"],
        nb=os.path.join(
            config["pipeline_path"],
            "notebooks/analyze_variant_counts.ipynb",
        ),
    output:
        nb="results/notebooks/analyze_variant_counts.ipynb",
        avg_counts_plot=config["variant_avg_counts_plot"],
        avg_counts_csv=config["variant_avg_counts_csv"],
    params:
        config["min_avg_counts"],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "analyze_variant_counts.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


rule check_adequate_variant_counts:
    """Check samples not specified for `exclude_after_counts` have adequate counts."""
    input:
        avg_counts_csv=rules.analyze_variant_counts.output.avg_counts_csv,
        avg_counts_plot=rules.analyze_variant_counts.output.avg_counts_plot,
        nb=rules.analyze_variant_counts.output.nb,
    output:
        # create flag file if all counts adequate
        passed=touch(os.path.join(config["variant_counts_dir"], "adequate_counts.flag")),
    params:
        min_avg_counts=config["min_avg_counts"],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "check_adequate_variant_counts.txt"),
    script:
        "scripts/check_adequate_variant_counts.py"


rule prob_escape:
    """Compute probabilities of escape for variants."""
    input:
        ancient(rules.check_adequate_variant_counts.output.passed),
        gene_sequence_codon=config["gene_sequence_codon"],
        codon_variants=config["codon_variants"],
        antibody_selections=config["antibody_selections"],
        site_numbering_map=config["site_numbering_map"],
        variant_counts=lambda wc: expand(
            rules.variant_counts.output.counts,
            library_sample=antibody_selection_group_samples[
                wc.antibody_selection_group
            ],
        ),
    output:
        prob_escape=os.path.join(
            config["prob_escape_dir"], "{antibody_selection_group}_prob_escape.csv"
        ),
        neut_standard_fracs=os.path.join(
            config["prob_escape_dir"],
            "{antibody_selection_group}_neut_standard_fracs.csv",
        ),
        neutralization=os.path.join(
            config["prob_escape_dir"], "{antibody_selection_group}_neutralization.csv"
        ),
    params:
        library_samples=lambda wc: antibody_selection_group_samples[
            wc.antibody_selection_group
        ],
        libraries=lambda wc: {
            libsamp: library_sample_to_library[libsamp]
            for libsamp in antibody_selection_group_samples[
                wc.antibody_selection_group
            ]
        },
        samples=lambda wc: {
            libsamp: library_sample_to_sample[libsamp]
            for libsamp in antibody_selection_group_samples[
                wc.antibody_selection_group
            ]
        },
        min_neut_standard_frac=config["prob_escape_min_neut_standard_frac"],
        min_neut_standard_count=config["prob_escape_min_neut_standard_count"],
        min_no_antibody_frac=config["prob_escape_min_no_antibody_frac"],
        min_no_antibody_counts=config["prob_escape_min_no_antibody_counts"],
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "prob_escape_{antibody_selection_group}.txt"),
    script:
        "scripts/prob_escape.py"


rule analyze_prob_escape:
    """Compute probabilities escape for variants."""
    input:
        prob_escape_files,
        config["antibody_selections"],
        nb=os.path.join(config["pipeline_path"], "notebooks/analyze_prob_escape.ipynb"),
    output:
        nb="results/notebooks/analyze_prob_escape.ipynb",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "analyze_prob_escape.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


rule fit_polyclonal:
    """Fit ``polyclonal`` models."""
    input:
        config["polyclonal_config"],
        prob_escape_csv=rules.prob_escape.output.prob_escape,
        nb=os.path.join(config["pipeline_path"], "notebooks/fit_polyclonal.ipynb"),
    output:
        pickle=os.path.join(
            config["polyclonal_dir"], "{antibody_selection_group}.pickle"
        ),
        nb="results/notebooks/fit_polyclonal_{antibody_selection_group}.ipynb",
    threads: config["fit_polyclonal_threads"]
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "fit_polyclonal_{antibody_selection_group}.txt"),
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p prob_escape_csv {input.prob_escape_csv} \
            -p pickle_file {output.pickle} \
            -p n_threads {threads} \
            &> {log}
        """
