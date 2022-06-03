"""Get variant counts from barcode counts."""


import sys

import alignparse.utils

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

barcode_runs = pd.read_csv(snakemake.input.barcode_runs)

variant_counts = pd.concat([
    pd.read_csv(variant_counts, na_filter=None).assign(library_sample=library_sample)
    for variant_counts, library_sample in
    zip(snakemake.input.variant_counts, snakemake.params.library_samples)
]).merge(barcode_runs[["library", "sample", "library_sample"]],
         on="library_sample", validate="many_to_one", how="left")
assert variant_counts.notnull().all().all()

variants.add_sample_counts_df(variant_counts)

antibody_selections = (
    pd.read_csv(snakemake.input.antibody_selections)
    .query("selection_group == @snakemake.wildcards.antibody_selection_group")
)

prob_escape, neut_standard_fracs, neutralization = variants.prob_escape(
    selections_df=antibody_selections,
    min_neut_standard_frac=snakemake.config["prob_escape_min_neut_standard_frac"],
    min_neut_standard_count=snakemake.config["prob_escape_min_neut_standard_count"],
)

# get the no-antibody count threshold and flag which prob_escape values meet it
threshold = (
    prob_escape.groupby(["library", "antibody_sample", "no-antibody_sample"], as_index=False)
    .aggregate(total_no_antibody_count=pd.NamedAgg("no-antibody_count", "sum"))
    .assign(
        no_antibody_count_threshold=lambda x: (
            x["total_no_antibody_count"] * snakemake.config["prob_escape_min_no_antibody_frac"]
        ).clip(lower=snakemake.config["prob_escape_min_no_antibody_counts"]).round().astype(int)
    )
)
prob_escape = prob_escape.merge(
    threshold, on=["library", "antibody_sample", "no-antibody_sample"], how="left", validate="many_to_one",
)
assert prob_escape["no_antibody_count_threshold"].notnull().all()

# renumber aa_substitutions from sequential to refernece numbering, also keep sequential
# numbering in column `aa_substitutions_sequential`. Also drop a few unneeded columns.
renumber = alignparse.utils.MutationRenumber(
    number_mapping=pd.read_csv(snakemake.input.site_numbering_map),
    old_num_col="sequential_site",
    new_num_col="reference_site",
    wt_nt_col=None,
)
prob_escape = (
    prob_escape.drop(
        columns=["codon_substitutions", "n_codon_substitutions"]
    )
    .rename(columns={"aa_substitutions": "aa_substitutions_sequential"})
    .assign(
        aa_substitutions=lambda x: (
            x["aa_substitutions_sequential"].apply(
                renumber.renumber_muts,
                allow_gaps=True,
                allow_stop=True,
            )
        ),
    )
    .merge(antibody_selections[["library", "antibody_sample", "no-antibody_sample", "antibody", "antibody_concentration"]], how="left", validate="many_to_one")
)

prob_escape.to_csv(snakemake.output.prob_escape, index=False, float_format="%.4f", na_rep="nan")
neut_standard_fracs.to_csv(snakemake.output.neut_standard_fracs, index=False, float_format="%.4g")
neutralization.to_csv(snakemake.output.neutralization, index=False, float_format="%.4g")
