"""Get functional scores from variant counts."""


import sys

import Bio.SeqIO

import alignparse.utils

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

variant_counts = pd.concat(
    [
        pd.read_csv(snakemake.input.preselection, na_filter=False).assign(
            library=snakemake.params.library,
            sample=snakemake.params.preselection_sample,
        ),
        pd.read_csv(snakemake.input.postselection, na_filter=False).assign(
            library=snakemake.params.library,
            sample=snakemake.params.postselection_sample,
        ),
    ]
)
assert variant_counts.notnull().all().all()

variants.add_sample_counts_df(variant_counts)

func_scores = variants.func_scores(
    preselection=snakemake.params.preselection_sample,
    pseudocount=snakemake.params.pseudocount,
    libraries=[snakemake.params.library],
)

# Renumber aa_substitutions in reference numbering, also keep sequential
# numbering in column `aa_substitutions_sequential`. Also drop a few unneeded columns.
# Also only keep the primary target.
renumber = alignparse.utils.MutationRenumber(
    number_mapping=pd.read_csv(snakemake.input.site_numbering_map),
    old_num_col="sequential_site",
    new_num_col="reference_site",
    wt_nt_col=None,
    allow_letter_suffixed_numbers=True,
)

func_scores = (
    func_scores.query("target == 'gene'")
    .drop(columns=["codon_substitutions", "n_codon_substitutions", "target"])
    .rename(columns={"aa_substitutions": "aa_substitutions_sequential"})
    .assign(
        aa_substitutions_reference=lambda x: (
            x["aa_substitutions_sequential"].apply(
                renumber.renumber_muts,
                allow_gaps=True,
                allow_stop=True,
            )
        ),
    )
)

func_scores.to_csv(
    snakemake.output.func_scores, index=False, float_format="%.4f", na_rep="nan"
)