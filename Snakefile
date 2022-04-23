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

rule align_parse_PacBio_ccs:
    """Align and parse PacBio CCS FASTQ file."""
    input:
        fastq=lambda wc: pacbio_runs.at[wc.pacbioRun, "fastq"],
        amplicon=config["pacbio_amplicon"],
        specs=config["pacbio_amplicon_specs"],
    output: outdir=directory(os.path.join(config["process_ccs_dir"], "{pacbioRun}"))
    conda: "environment.yml"
    script: "scripts/align_parse_PacBio_ccs.py"

rule analyze_pacbio_ccs:
    """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
    input:
        expand(rules.align_parse_PacBio_ccs.output.outdir, pacbioRun=pacbio_runs.index),
        config["pacbio_amplicon"],
        config["pacbio_amplicon_specs"],
        nb=os.path.join(config["pipeline_path"], "notebooks/analyze_pacbio_ccs.ipynb"),
    output:
        config["aligned_ccs_file"],
        nb="results/notebooks/analyze_pacbio_ccs.ipynb"
    conda: "environment.yml"
    shell: "papermill {input.nb} {output.nb}"
