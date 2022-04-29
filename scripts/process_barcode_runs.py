"""Process barcode runs to add samples and ensure runs are unique."""


import os
import sys

import pandas as pd


input_csv = snakemake.input.csv
output_csv = snakemake.output.csv
pacbio_libraries = snakemake.params.pacbio_libraries


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")


def process_sample(row):
    if row["library"] not in pacbio_libraries:
        raise ValueError(f"library {row['library']} not in {pacbio_libraries=}")
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


print(f"Reading barcode runs from {input_csv}")
barcode_runs = pd.read_csv(snakemake.input.csv).assign(
    sample=lambda x: x.apply(process_sample, axis=1),
    library_sample=lambda x: x["library"] + "_" + x["sample"],
)

# ensure samples are unique
dups = (
    barcode_runs.groupby("library_sample")
    .aggregate(n=pd.NamedAgg("library", "count"))
    .query("n > 1")
)
if len(dups):
    raise ValueError(f"Found some duplicated samples:\n{dups}")

# ensure FASTQs exist and are unique
fastqs = (
    barcode_runs.assign(
        fastq=lambda x: x["fastq_R1"].map(lambda fs: [f.strip() for f in fs.split(";")])
    )
    .explode("fastq")
    .assign(
        found_file=lambda x: x["fastq"].map(os.path.isfile),
        n_occurrences=lambda x: x.groupby("fastq")["library_sample"].transform("count"),
    )[["library_sample", "fastq", "found_file", "n_occurrences"]]
)
if not fastqs["found_file"].all():
    raise ValueError(f"Failed to find some fastqs:\n{fastqs.query('not found_file')}")
if any(fastqs["n_occurrences"] != 1):
    raise ValueError(f"FASTQs repeated:\n{fastqs.query('n_occurrences != 1')}")

print(f"Writing barcode runs with samples to {output_csv}")
barcode_runs.to_csv(output_csv, index=False)
