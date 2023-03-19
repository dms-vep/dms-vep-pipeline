"""Get probabilities of escape from variant counts."""


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

library = set(snakemake.params.libraries.values())
assert len(library) == 1
library = list(library)[0]

antibody = set(snakemake.params.antibodies)
assert len(antibody) == 1
antibody = list(antibody)[0]

variant_counts = pd.concat(
    [
        pd.read_csv(variant_counts, na_filter=None).assign(
            library_sample=library_sample
        )
        for variant_counts, library_sample in zip(
            snakemake.input.variant_counts, snakemake.params.library_samples
        )
    ]
).assign(
    library=library,
    sample=lambda x: x["library_sample"].map(snakemake.params.samples),
)
assert variant_counts.notnull().all().all()

variants.add_sample_counts_df(variant_counts)

selections_df = pd.DataFrame(
    {
        "library": library,
        "antibody_sample": snakemake.params.antibody_samples,
        "no-antibody_sample": snakemake.params.no_antibody_samples,
        "antibody": antibody,
        "antibody_concentration": snakemake.params.antibody_concentrations,
    }
)

neut_standard = sorted(set(snakemake.params.neut_standard))
if len(neut_standard) != 1:
    raise ValueError(f"{neut_standard=} not unique for\n{selections_df=}")
neut_standard = neut_standard[0]

prob_escape, neut_standard_fracs, neutralization = variants.prob_escape(
    selections_df=selections_df,
    neut_standard_target=neut_standard,
    min_neut_standard_frac=snakemake.params.min_neut_standard_frac,
    min_neut_standard_count=snakemake.params.min_neut_standard_count,
)

# get the no-antibody count threshold and flag which prob_escape values meet it
threshold = (
    prob_escape.groupby(
        ["library", "antibody_sample", "no-antibody_sample"], as_index=False
    )
    .aggregate(total_no_antibody_count=pd.NamedAgg("no-antibody_count", "sum"))
    .assign(
        no_antibody_count_threshold=lambda x: (
            x["total_no_antibody_count"] * snakemake.params.min_no_antibody_frac
        )
        .clip(lower=snakemake.params.min_no_antibody_counts)
        .round()
        .astype(int)
    )
)
prob_escape = prob_escape.merge(
    threshold,
    on=["library", "antibody_sample", "no-antibody_sample"],
    how="left",
    validate="many_to_one",
)
assert prob_escape["no_antibody_count_threshold"].notnull().all()

# renumber aa_substitutions in reference numbering, also keep sequential
# numbering in column `aa_substitutions_sequential`. Also drop a few unneeded columns.
renumber = alignparse.utils.MutationRenumber(
    number_mapping=pd.read_csv(snakemake.input.site_numbering_map),
    old_num_col="sequential_site",
    new_num_col="reference_site",
    wt_nt_col=None,
    allow_letter_suffixed_numbers=True,
)
prob_escape = (
    prob_escape.drop(columns=["codon_substitutions", "n_codon_substitutions"])
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
    .merge(selections_df, how="left", validate="many_to_one")
)

prob_escape.to_csv(
    snakemake.output.prob_escape, index=False, float_format="%.4f", na_rep="nan"
)
neut_standard_fracs.to_csv(
    snakemake.output.neut_standard_fracs, index=False, float_format="%.4g"
)
neutralization.to_csv(snakemake.output.neutralization, index=False, float_format="%.4g")
