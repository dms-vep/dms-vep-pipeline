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

    def process_exclude(s):
        """Process a string specifing whether to exclude."""
        if pd.isnull(s) or s.lower() in {"no", "false"}:
            return "no"
        elif s.lower() in {"yes", "true"}:
            return "yes"

    barcode_runs = pd.read_csv(barcode_runs_csv).assign(
        sample=lambda x: x.apply(process_sample, axis=1),
        library_sample=lambda x: x["library"] + "_" + x["sample"],
        fastq_R1=lambda x: x["fastq_R1"].map(
            lambda fs: [f.strip() for f in fs.split(";")]
        ),
        exclude_after_counts=lambda x: x["exclude_after_counts"].map(process_exclude),
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


def get_antibody_selections(
    barcode_runs,
    pair_on=("library", "virus_batch", "date", "replicate"),
):
    """Data frame of antibody selections from data frame of barcode runs.

    Pairs each antibody selection with the corresponding no-antibody control
    that is the same in all of the properties specified in `pair_on`.

    """
    barcode_runs = barcode_runs.query("exclude_after_counts == 'no'")

    antibodies = barcode_runs.query("sample_type == 'antibody'").rename(
        columns={"sample": "antibody_sample"}
    )[["antibody_sample", "antibody_concentration", *pair_on]]
    assert len(antibodies) == len(antibodies.drop_duplicates())

    controls = barcode_runs.query("sample_type == 'no-antibody_control'").rename(
        columns={"sample": "no-antibody_sample"}
    )[["no-antibody_sample", *pair_on]]
    assert len(controls) == len(controls.drop_duplicates())

    antibody_selections = (
        antibodies.merge(controls, how="left", validate="many_to_one", on=pair_on)
        .merge(
            barcode_runs[["library", "sample", "library_sample"]].rename(
                columns={
                    "sample": "antibody_sample",
                    "library_sample": "antibody_library_sample",
                },
            ),
            how="left",
            validate="many_to_one",
            on=["library", "antibody_sample"],
        )
        .merge(
            barcode_runs[["library", "sample", "library_sample"]].rename(
                columns={
                    "sample": "no-antibody_sample",
                    "library_sample": "no-antibody_library_sample",
                },
            ),
            how="left",
            validate="many_to_one",
            on=["library", "no-antibody_sample"],
        )
    )

    if antibody_selections.isnull().any().any():
        raise ValueError(
            "null antibody selections:\n"
            + str(antibody_selections[antibody_selections.isnull().any(axis=1)])
        )

    assert (
        len(antibody_selections)
        == len(antibody_selections.groupby(["library", "antibody_sample"]))
        == antibody_selections["antibody_library_sample"].nunique()
    )

    return antibody_selections
