# Input data
This subdirectory contains input data used by the pipeline.

## PacBio full-length variant sequencing to link barcodes

[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file having features to parse with [alignparse](https://jbloomlab.github.io/alignparse/). Must have *gene* (the gene of interest) and *barcode* features.

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/).

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants.

## Illumina barcode sequencing

[barcode_runs.csv](barcode_runs.csv): List of Illumina barcode runs.

## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential 1, 2, ... numbering of the gene to a "reference" numbering scheme that represents the standard naming of sites for this gene.

## Neutralization standard barcodes
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) barcodes for the neutralization standard.

## Configuration for `polyclonal` fitting
[polyclonal_config.yaml](polyclonal_config.yaml) specifies how the analysis with [polyclonal](https://jbloomlab.github.io/polyclonal/) is done.
For each antibody listed in [barcode_runs.csv](barcode_runs.csv), specify:

 - *max_epitopes*: the maximum number of epitopes to test. The fitting keeps testing more epitopes up to this max until additional epitopes don't improve fitting anymore.
 - *n_bootstrap_samples*: number of bootstrap samples to use.
 - *reg_escape_weight*: regularization weight for mutation-escape values.
 - *reg_spread_weight*: regularization weight for spread of escape values at each site.
 - *reg_activity_weight*: regularization weight for epitope activities.
 - *times_seen*: the `times_seen` value used for plotting the results (number of variants a mutation must be found in).
 - *min_epitope_activity_to_include*: keep adding epitopes until activity <= this.
