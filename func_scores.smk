"""``snakemake`` file for functional scores steps of pipeline."""

if "prebuilt_muteffects" in config and config["prebuilt_muteffects"]:
    # Use pre-built mutational effects
    have_muteffects = True
    func_selections = pd.DataFrame()  # make empty data frame as no func selections
    muteffects_plots = {}

    rule get_muteffects:
        """Get prebuilt mutational effects."""
        output:
            muteffects=config["muteffects_observed"],
        params:
            url=config["prebuilt_muteffects"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "get_muteffects.txt"),
        shell:
            "curl -o {output.muteffects} {params.url} &> {log}"

else:
    # Compute functional scores and mutational effects de novo
    func_selections = get_functional_selections(barcode_runs)
    os.makedirs(os.path.dirname(config["functional_selections"]), exist_ok=True)
    to_csv_if_changed(func_selections, config["functional_selections"], index=False)
    func_selections_dict = func_selections.set_index("selection_name").to_dict(
        orient="index"
    )

    func_score_files = [
        os.path.join(config["func_score_dir"], f"{func_selection}_func_scores.csv")
        for func_selection in func_selections["selection_name"]
    ]

    if len(func_selections):
        muteffects_plots = {
            f"muteffects_{pheno}_heatmap": os.path.splitext(
                config[f"muteffects_{pheno}"]
            )[0]
            + "_heatmap.html"
            for pheno in ["observed", "latent"]
        }
        have_muteffects = True
    else:
        muteffects_plots = {}
        have_muteffects = False

    muteffects_files = [
        os.path.join(
            config["globalepistasis_dir"],
            f"{func_selection}_muteffects_{phenotype}.csv",
        )
        for func_selection in func_selections["selection_name"]
        for phenotype in ["latent", "observed"]
    ]

    rule func_scores:
        """Compute functional scores for variants."""
        input:
            gene_sequence_codon=config["gene_sequence_codon"],
            codon_variants=config["codon_variants"],
            site_numbering_map=config["site_numbering_map"],
            preselection=lambda wc: os.path.join(
                config["variant_counts_dir"],
                func_selections_dict[wc.func_selection]["preselection_library_sample"]
                + ".csv",
            ),
            postselection=lambda wc: os.path.join(
                config["variant_counts_dir"],
                func_selections_dict[wc.func_selection]["postselection_library_sample"]
                + ".csv",
            ),
        output:
            func_scores=os.path.join(
                config["func_score_dir"],
                "{func_selection}_func_scores.csv",
            ),
        params:
            library=lambda wc: func_selections_dict[wc.func_selection]["library"],
            preselection_sample=lambda wc: func_selections_dict[wc.func_selection][
                "preselection_sample"
            ],
            postselection_sample=lambda wc: func_selections_dict[wc.func_selection][
                "postselection_sample"
            ],
            pseudocount=config["func_scores_pseudocount"],
            min_wt_count=config["func_scores_min_wt_count"],
            min_wt_frac=config["func_scores_min_wt_frac"],
            min_preselection_counts=config["func_scores_min_preselection_counts"],
            min_preselection_frac=config["func_scores_min_preselection_frac"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "func_scores_{func_selection}.txt"),
        script:
            "scripts/func_scores.py"

    rule analyze_func_scores:
        """Analyze the functional scores."""
        input:
            func_score_files,
            config["functional_selections"],
            config["mutation_design_classification"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/analyze_func_scores.ipynb"
            ),
        output:
            # only make a notebook output for docs if there are functional selections
            **(
                {"nb": "results/notebooks/analyze_func_scores.ipynb"}
                if len(func_selections)
                else {}
            ),
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "analyze_func_scores.txt"),
        shell:
            "papermill {input.nb} {output.nb} &> {log}"

    rule fit_globalepistasis:
        """Fit global epistasis models to variant functional scores to get muteffects."""
        input:
            func_scores_csv=rules.func_scores.output.func_scores,
            site_numbering_map=config["site_numbering_map"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/fit_globalepistasis.ipynb"
            ),
        output:
            pickle=os.path.join(
                config["globalepistasis_dir"], "{func_selection}.pickle"
            ),
            muteffects_latent=os.path.join(
                config["globalepistasis_dir"],
                "{func_selection}_muteffects_latent.csv",
            ),
            muteffects_observed=os.path.join(
                config["globalepistasis_dir"],
                "{func_selection}_muteffects_observed.csv",
            ),
            nb="results/notebooks/fit_globalepistasis_{func_selection}.ipynb",
        params:
            func_scores_floor=(
                config["func_scores_floor"] if "func_scores_floor" in config else None
            ),
            plot_kwargs_yaml=yaml.dump(
                {"plot_kwargs": config["muteffects_plot_kwargs"]}
            ),
            likelihood=(
                config["epistasis_model_likelihood"]
                if "epistasis_model_likelihood" in config
                else "Gaussian"
            ),
            ftol=(
                config["epistasis_model_ftol"]
                if "epistasis_model_ftol" in config
                else 1e-7
            ),
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "fit_globalepistasis_{func_selection}.txt"),
        shell:
            """
            papermill {input.nb} {output.nb} \
                -p func_scores_csv {input.func_scores_csv} \
                -p sitenumbering_map_csv {input.site_numbering_map} \
                -p pickle_file {output.pickle} \
                -p muteffects_latent_csv {output.muteffects_latent} \
                -p muteffects_observed_csv {output.muteffects_observed} \
                -p func_scores_floor {params.func_scores_floor} \
                -p likelihood {params.likelihood} \
                -p ftol {params.ftol} \
                -y "{params.plot_kwargs_yaml}" \
                &> {log}
            """

    rule avg_muteffects:
        """Average the mutation effects on viral entry across replicates and libraries."""
        input:
            config["functional_selections"],
            expand(
                rules.fit_globalepistasis.output.muteffects_latent,
                func_selection=func_selections["selection_name"],
            ),
            expand(
                rules.fit_globalepistasis.output.muteffects_observed,
                func_selection=func_selections["selection_name"],
            ),
            nb=os.path.join(config["pipeline_path"], "notebooks/avg_muteffects.ipynb"),
        output:
            config["muteffects_observed"],
            config["muteffects_latent"],
            os.path.splitext(config["muteffects_observed"])[0]
            + "_heatmap_unformatted.html",
            os.path.splitext(config["muteffects_latent"])[0]
            + "_heatmap_unformatted.html",
            # only make a notebook output for docs if there are functional selections
            **(
                {"nb": "results/notebooks/avg_muteffects.ipynb"}
                if len(func_selections)
                else {}
            ),
        params:
            config["muteffects_plot_kwargs"],
            config["muteffects_avg_method"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "avg_muteffects.txt"),
        shell:
            "papermill {input.nb} {output.nb} &> {log}"

    rule format_muteffects_plot:
        """Format muteffects plot."""
        input:
            chart=(
                "results/muteffects_functional/muteffects_{pheno}_heatmap_unformatted.html"
            ),
            md=(
                config["muteffects_legend"]
                if "muteffects_legend" in config and config["muteffects_legend"]
                else os.path.join(
                    config["pipeline_path"], "plot_legends/muteffects_legend.md"
                )
            ),
            pyscript=os.path.join(
                config["pipeline_path"], "scripts/format_altair_html.py"
            ),
        output:
            chart="results/muteffects_functional/muteffects_{pheno}_heatmap.html",
        params:
            format_plot=int(
                ("format_muteffects_plots" not in config)
                or config["format_muteffects_plots"]
            ),
            site=lambda _, output: os.path.join(
                github_pages_url,
                os.path.basename(output.chart),
            ),
            title=f"mutation effects for {config['github_repo']}",
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "format_muteffects_plot_{pheno}.txt"),
        shell:
            """
            if [[ {params.format_plot} -eq 1 ]]; then
                python {input.pyscript} \
                    --chart {input.chart} \
                    --markdown {input.md} \
                    --site "{params.site}" \
                    --title "{params.title}" \
                    --description "Mutational effects on {wildcards.pheno} phenotype" \
                    --output {output} \
                    &> {log}
            else
                cp {input.chart} {output.chart}
            fi
            """
