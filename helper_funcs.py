"""Helper functions for the ``snakemake`` pipeline."""


import os

import pandas as pd


def pacbio_runs_from_config(pacbio_runs_csv):
    """Data frame of PacBio runs from input CSV."""
    pacbio_runs = (
        pd.read_csv(pacbio_runs_csv, dtype=str)
        .assign(pacbioRun=lambda x: x["library"] + "_" + x["run"].astype(str))
        .set_index("pacbioRun")
    )
    assert pacbio_runs.index.nunique() == len(pacbio_runs)
    return pacbio_runs


def barcode_runs_from_config(barcode_runs_csv, valid_libraries):
    """Data frame of barcode runs from input CSV."""

    def process_sample(row):
        """Internal function that processes rows in data frame."""
        if row["library"] not in valid_libraries:
            raise ValueError(f"library {row['library']} not in {valid_libraries=}")
        label_cols = ["date", "virus_batch", "sample_type"]
        if row["sample_type"] == "antibody":
            label_cols += ["antibody", "antibody_concentration"]
        elif row["sample_type"] not in {"VSVG_control", "no-antibody_control"}:
            raise ValueError(f"invalid `sample_type` {row['sample_type']}")
        label_cols.append("replicate")
        sample = []
        for col in label_cols:
            if pd.isnull(row[col]):
                raise ValueError(f"null {col} in {row}")
            sample.append(row[col])
        return "_".join(map(str, sample))

    barcode_runs = pd.read_csv(barcode_runs_csv).assign(
        sample=lambda x: x.apply(process_sample, axis=1),
        library_sample=lambda x: x["library"] + "_" + x["sample"],
        fastq_R1=lambda x: x["fastq_R1"].map(
            lambda fs: [f.strip() for f in fs.split(";")]
        ),
    )

    # check no duplicate samples
    dups = (
        barcode_runs.groupby("library_sample")
        .aggregate(n=pd.NamedAgg("library", "count"))
        .query("n > 1")
    )
    if len(dups):
        raise ValueError(f"Found some duplicated samples:\n{dups}")

    # ensure FASTQs exist and are unique
    fastqs = barcode_runs.explode("fastq_R1").assign(
        found_file=lambda x: x["fastq_R1"].map(os.path.isfile),
        n_occurrences=lambda x: x.groupby("fastq_R1")["library_sample"].transform(
            "count"
        ),
    )[["library_sample", "fastq_R1", "found_file", "n_occurrences"]]
    if not fastqs["found_file"].all():
        raise ValueError(
            f"Failed to find some fastqs:\n{fastqs.query('not found_file')}"
        )
    dup_fastqs = fastqs.query("n_occurrences != 1")
    if any(fastqs["n_occurrences"] != 1):
        pd.set_option("display.max_colwidth", None)
        raise ValueError(
            f"FASTQs repeated:\n{dup_fastqs[['fastq_R1', 'n_occurrences']]}"
        )

    return barcode_runs
