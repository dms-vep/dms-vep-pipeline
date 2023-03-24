"""``snakemake`` file that includes the pipeline code.

It doesn't include the `configfile`, which is designed to be contained
in the upstream file that includes this one.

"""

# Imports ---------------------------------------------------------------------
import glob
import os
import re

import yaml


include: "funcs.smk"  # import functions


# Check to make sure Github repo information correct in config ----------------
if ("no_github_check" not in config) or (not config["no_github_check"]):
    git_remote_res = subprocess.run(
        ["git", "remote", "-v"],
        capture_output=True,
        text=True,
    )
    git_remote_regex = re.match(
        "origin\t(https://|git@)github.com(/|:)(?P<user>[\-\w]+)/(?P<repo>[\-\w\.]+)(?: |\.git)",
        git_remote_res.stdout,
    )
    if not git_remote_regex:
        raise ValueError(f"cannot match git repo from\n{git_remote_res}")
    for attr in ["repo", "user"]:
        regex_attr = git_remote_regex.group(attr)
        if attr == "repo" and regex_attr.endswith(".git"):
            regex_attr = regex_attr[:-4]
        config_attr = config[f"github_{attr}"]
        if regex_attr != config_attr:
            raise ValueError(
                f"github_{attr} in `config.yaml` does not match actual git remote:\n"
                f"{regex_attr} versus {config_attr}"
            )

github_pages_url = f"https://{config['github_user']}.github.io/{config['github_repo']}"


# Global variables and processing before pipeline -----------------------------

# Data frames for PacBio runs, Illumina barcode runs, antibody selections, etc.
# Some of these are written to CSV files, but only if they have changed.

barcode_runs = barcode_runs_from_config(config["barcode_runs"])
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

antibody_escape_files = [
    os.path.join(config["escape_dir"], f"{antibody}_{suffix}")
    for antibody in antibody_selections["antibody"].unique()
    for suffix in ["avg.csv", "rep.csv"]
]

antibody_escape_plots = [
    os.path.join(config["escape_dir"], f"{antibody}_escape_plot.html")
    for antibody in antibody_selections["antibody"].unique()
]


# Rules ---------------------------------------------------------------------

include: "build_variants.smk"  # build variants with included rules
include: "func_scores.smk"  # get functional scores with included rules

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
        # only make a notebook output for docs if there are barcode runs
        **(
            {"nb": "results/notebooks/analyze_variant_counts.ipynb"}
            if len(barcode_runs)
            else {}
        ),
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
        rules.analyze_variant_counts.output.avg_counts_csv if len(barcode_runs) else [],
    output:
        # create flag file if all counts adequate
        passed=touch(os.path.join(config["variant_counts_dir"], "adequate_counts.flag")),
    params:
        min_avg_counts=config["min_avg_counts"],
        barcode_runs_exist=(len(barcode_runs) > 0),
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
        antibody_samples=lambda wc: tuple(
            antibody_selections.query(
                "selection_group == @wc.antibody_selection_group"
            )["antibody_sample"]
        ),
        no_antibody_samples=lambda wc: tuple(
            antibody_selections.query(
                "selection_group == @wc.antibody_selection_group"
            )["no-antibody_sample"]
        ),
        antibodies=lambda wc: tuple(
            antibody_selections.query(
                "selection_group == @wc.antibody_selection_group"
            )["antibody"]
        ),
        antibody_concentrations=lambda wc: tuple(
            antibody_selections.query(
                "selection_group == @wc.antibody_selection_group"
            )["antibody_concentration"]
        ),
        neut_standard=lambda wc: tuple(
            antibody_selections.query(
                "selection_group == @wc.antibody_selection_group"
            )["neut_standard_name"]
        ),
        min_neut_standard_frac=config["prob_escape_min_neut_standard_frac"],
        min_neut_standard_count=config["prob_escape_min_neut_standard_count"],
        min_no_antibody_frac=config["prob_escape_min_no_antibody_frac"],
        min_no_antibody_counts=config["prob_escape_min_no_antibody_counts"],
        min_antibody_frac=(
            config["prob_escape_min_antibody_frac"]
            if "prob_escape_min_antibody_frac" in config
            else None
        ),
        min_antibody_counts=(
            config["prob_escape_min_antibody_counts"]
            if "prob_escape_min_antibody_counts" in config
            else None
        ),
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
        # only make a notebook output for docs if there are antibody selections
        **(
            {"nb": "results/notebooks/analyze_prob_escape.ipynb"}
            if len(antibody_selections)
            else {}
        ),
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
        config["site_numbering_map"],
        **(
            {"spatial_distances": config["spatial_distances"]}
            if "spatial_distances" in config
            else {}
        ),
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


rule avg_antibody_escape:
    """Average escape for an antibody or serum."""
    input:
        **({"muteffects": config["muteffects_observed"]} if have_muteffects else {}),
        site_numbering_map=config["site_numbering_map"],
        polyclonal_config=config["polyclonal_config"],
        selection_group_pickles=lambda wc: expand(
            rules.fit_polyclonal.output.pickle,
            antibody_selection_group=(
                antibody_selections.query("antibody == @wc.antibody")[
                    "selection_group"
                ].unique()
            ),
        ),
        nb=os.path.join(config["pipeline_path"], "notebooks/avg_antibody_escape.ipynb"),
    output:
        avg_pickle=os.path.join(config["escape_dir"], "{antibody}.pickle"),
        avg_escape=os.path.join(config["escape_dir"], "{antibody}_avg.csv"),
        rep_escape=os.path.join(config["escape_dir"], "{antibody}_rep.csv"),
        escape_plot=os.path.join(
            config["escape_dir"],
            "{antibody}_escape_plot_unformatted.html",
        ),
        nb="results/notebooks/avg_antibody_escape_{antibody}.ipynb",
    params:
        escape_avg_method=config["escape_avg_method"],
        selection_groups_yaml=lambda wc: yaml.dump(
            {
                "selection_groups_dict": (
        antibody_selections.query("antibody == @wc.antibody")[
            [
                "library",
                "virus_batch",
                "date",
                "replicate",
                "selection_group",
            ]
        ]
        .drop_duplicates()
        .assign(
            pickle_file=lambda x: (
                config["polyclonal_dir"]
        + "/"
                            + x["selection_group"]
                            + ".pickle"
                        )
                    )
                    .set_index("selection_group")
                    .to_dict(orient="index")
                )
            }
        ),
        muteffects=lambda _, input: input.muteffects if have_muteffects else "none",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "avg_antibody_escape_{antibody}.txt"),
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p antibody {wildcards.antibody} \
            -p escape_avg_method {params.escape_avg_method} \
            -p polyclonal_config {input.polyclonal_config} \
            -p muteffects_csv {params.muteffects} \
            -p site_numbering_map {input.site_numbering_map} \
            -p avg_pickle {output.avg_pickle} \
            -p avg_escape {output.avg_escape} \
            -p rep_escape {output.rep_escape} \
            -p escape_plot {output.escape_plot} \
            -y "{params.selection_groups_yaml}" \
            &> {log}
        """


rule format_antibody_escape_plot:
    """Add formatting to antibody escape plots."""
    input:
        chart=rules.avg_antibody_escape.output.escape_plot,
        md=(
            config["antibody_escape_legend"]
            if "antibody_escape_legend" in config and config["antibody_escape_legend"]
            else os.path.join(
                config["pipeline_path"], "plot_legends/antibody_escape_legend.md"
            )
        ),
        pyscript=os.path.join(config["pipeline_path"], "scripts/format_altair_html.py"),
    output:
        chart=os.path.join(
            config["escape_dir"],
            "{antibody}_escape_plot.html",
        ),
    params:
        format_plot=int(
            ("format_antibody_escape_plots" not in config)
            or config["format_antibody_escape_plots"]
        ),
        site=lambda _, output: os.path.join(
            github_pages_url,
            os.path.basename(output.chart),
        ),
        title=lambda wc: f"{wc.antibody} for {config['github_repo']}",
    conda:
        "environment.yml"
    log:
        os.path.join(config["logdir"], "format_antibody_escape_plot_{antibody}.txt"),
    shell:
        """
        if [[ {params.format_plot} -eq 1 ]]; then
            python {input.pyscript} \
                --chart {input.chart} \
                --markdown {input.md} \
                --site "{params.site}" \
                --title "{params.title}" \
                --description "Interactive plot of antibody escape" \
                --output {output} \
                &> {log}
        else
            cp {input.chart} {output.chart}
        fi
        """
