"""Get variant counts from barcode counts."""


import sys

import Bio.SeqIO

import dms_variants.codonvarianttable

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

geneseq = str(Bio.SeqIO.read(snakemake.input.gene_sequence_codon, "fasta").seq)

variants = dms_variants.codonvarianttable.CodonVariantTable(
    barcode_variant_file=snakemake.input.codon_variants,
    geneseq=geneseq,
    allowgaps=True,
    substitutions_are_codon=True,
    primary_target="gene",
    substitutions_col="codon_substitutions",
)

variants.addSampleCounts(
    library=snakemake.params.library,
    sample=snakemake.params.sample,
    barcodecounts=pd.read_csv(snakemake.input.barcode_counts),
)

(
    variants.variant_count_df[
        [
            "barcode",
            "count",
            "codon_substitutions",
            "aa_substitutions",
            "variant_call_support",
        ]
    ]
    .sort_values(["count", "barcode"], ascending=[False, True])
    .to_csv(snakemake.output.counts, index=False)
)
