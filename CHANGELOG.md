# Change log

### version 2.2.0
- Add `epitope_colors` option in `polyclonal_config` to specify epitope colors.

#### version 2.1.1
- Update to `polyclonal` 5.1 and `altair` 5.0.0rc1

### version 2.1.0
- Added `environment_align_parse_PacBio_ccs.yml` environment specific to the `align_parse_PacBio_ccs` rule so that long-running rule isn't re-run every time main environment in `environment.yml` is updated.
- Update to `polyclonal` version 5.0. This will again change results of fitting models.

#### version 2.0.1
- Update `polyclonal` to 4.1.

## version 2.0.0
- Update `polyclonal` to version 4.0. This is a new major version of `dms-vep-pipeline` as it will change the results of running model fits if default parameters are used (they can always be adjusted by setting the keyword arguments to `Polyclonal.fit` via `fit_kwargs` in `polyclonal_config`.) See [here](https://github.com/jbloomlab/polyclonal/pull/151#issuecomment-1440892895) for more details.
- Update notebooks with `polyclonal` fitting results to show `curves_plot` rather than `activity_wt_barplot`.

#### version 1.8.1
- Fix bug in `avg_antibody_escape.ipynb` that was making it always take the model with fewest epitopes. Now takes model with most that can harmonize as intended.

### version 1.8.0
- Upgrade to `polyclonal` version 3.4. This version has improved plotting including showing individual libraries on tooltips, better handling NaN values in heatmap tooltips, and hiding rather than filtering mutations based on functional effect.
- Also updated versions of other software in `conda` environment, although updates other than `polyclonal` should not change output.

### version 1.7.0
- add option to use pre-built mutational effects.
 - Enables use of functional scores from another repo, do not calculate here. Useful for small repos that for instance might add another antibody to an already extensively studied library.
 - Add `prebuilt_muteffects` option in `config.yaml`.
 - Move functional score computation from `pipeline.smk` to `func_scores.smk`

#### version 1.6.1
- Fix bug in using pre-built variants when no `pacbio_runs` file provided.

### version 1.6.0
- add option to use pre-built variants
 - Enables use of pre-built codon-variant table rather than having to rebuild every time.
 - Add `prebuilt_variants` and `prebuilt_geneseq` options
 - Move variant building from `pipeline.smk` to `build_variants.smk`

#### version 1.5.2
- Update `polyclonal` to 3.3

#### version 1.5.1
- Update `polyclonal` to 3.2

### version 1.5.0
- Better harmonization of epitopes in `polyclonal` fits with multiple epitopes:
 - The `fit_polyclonal` rule now pickles a list of all models (with one, two, etc epitopes) rather than the single model with most epitopes that fits stopping criteria.
 - The `avg_antibody_escape` rule now tries to harmonize the models with the most epitopes, and then keeps decreasing number of epitopes until they can be harmonized.

#### version 1.4.3
- Upgrade to `biopython` 1.80.

#### version 1.4.2
- Plot of functional scores in `analyze_func_scores` separates variants by whether they have only intended or also unintended nonsynonymous mutations.

#### version 1.4.1
- Update `upsetplot` and `matplotlib` in `conda` environment.

### version 1.4
- Add `spatial_distances` to `config.yaml` to enable spatial regularization by `polyclonal`.
- Changes to `conda` environment in `environment.yml`:
 - Upgrade to `polyclonal` version 3.1
 - Upgrade to `alignparse` version 0.6.0
 - Upgrade `mafft` to 7.508
 - Upgrade `snakemake` to 7.18.2
 - No longer pin `tabulate` version (see #90)
- Fix check for repeated FASTQs even when *notes* column in `barcode_runs.csv` is all *null*.

#### version 1.3.1
- A FASTQ can be repeated in `barcode_runs.csv` if the *notes* column has the word "repeated" for each repeated row.

### version 1.3
- Upgraded `polyclonal` to version 2.5.

### version 1.2
- Upgrade `python` from 3.9 to 3.10
- Improve text formatting in plots with text added by `format_altair_html.py`.
- Regex check for correct repo in snakemake pipeline checks ssh repo origins.
- Antibody escape analyses fixed to work even when no functional effects data.

### version 1.1
- Change handling of formatted antibody-escape and functional effect plots. Format the plots **by default** unless config specifies not to, and do not suffix the formatted plots.

#### version 1.01
- Update `polyclonal` to version 2.4 to fix problem with wildtype in heatmaps.
- Fix `tabulate` to <0.9 in `conda` environment due to this [issue](https://github.com/snakemake/snakemake/issues/1891)

## version 1.0
Initial commit
