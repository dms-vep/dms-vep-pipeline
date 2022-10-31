# Change log

### version 1.dev2
- Improve text formatting in plots with text added by `format_altair_html.py`.
- Regex check for correct repo in snakemake pipeline checks ssh repo origins.

### version 1.1
- Change handling of formatted antibody-escape and functional effect plots. Format the plots **by default** unless config specifies not to, and do not suffix the formatted plots.

#### version 1.01
- Update `polyclonal` to version 2.4 to fix problem with wildtype in heatmaps.
- Fix `tabulate` to <0.9 in `conda` environment due to this [issue](https://github.com/snakemake/snakemake/issues/1891)

## version 1.0
Initial commit
