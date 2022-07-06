"""Check that samples have adequate variant counts, raise informative error if not."""


import sys

import pandas as pd

orig_sys_stderr = sys.stderr
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

if not snakemake.params.barcode_runs_exist:
    print("No barcode runs")
    sys.exit(0)

avg_counts_csv = snakemake.input[0]
avg_counts = pd.read_csv(avg_counts_csv)

too_few = (
    avg_counts.query("(not adequate_counts) & (exclude_after_counts != 'yes')")
    .drop(columns=["adequate_counts", "min_avg_counts"])
    .reset_index(drop=True)
)

if len(too_few):
    err_msg = (
        f"\n{len(too_few)} samples have fewer than {snakemake.params.min_avg_counts} "
        "average counts per variant but do not have `exclude_after_counts` set to "
        "`yes` in `barcode_runs`. You need to exclude these samples after counts.\n\n"
        "The average counts per variant for each sample are tallied in {avg_counts_csv}"
        f"\n\nHere are the problematic samples:\n{str(too_few)}\n"
    )
    orig_sys_stderr.write(err_msg)
    orig_sys_stderr.flush()
    raise ValueError(err_msg)
else:
    print("All samples have adequate counts or have `exclude_after_counts` specified.")
