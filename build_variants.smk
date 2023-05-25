"""``snakemake`` file that includes the pipeline code for building variants."""


if "prebuilt_variants" in config and config["prebuilt_variants"]:
    if not (("prebuilt_geneseq" in config) and config["prebuilt_geneseq"]):
        raise ValueError(
            "specify both or neither of `prebuilt_geneseq` and `prebuilt_variants`"
        )

    rule get_prebuilt_variants:
        """Get pre-built variants and gene sequence."""
        output:
            variants=config["codon_variants"],
            geneseq=config["gene_sequence_codon"],
        params:
            variants_url=config["prebuilt_variants"],
            geneseq_url=config["prebuilt_geneseq"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "get_prebuilt_variants.txt"),
        shell:
            """
            curl -o {output.variants} {params.variants_url} &> {log}
            curl -o {output.geneseq} {params.geneseq_url} &> {log}
            """

    rule translate_geneseq:
        """Translate gene sequence into protein sequence."""
        input:
            gene=rules.get_prebuilt_variants.output.geneseq,
        output:
            prot=config["gene_sequence_protein"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "translate_geneseq.txt"),
        script:
            "scripts/translate_geneseq.py"

else:
    if ("prebuilt_geneseq" in config) and config["prebuilt_geneseq"]:
        raise ValueError(
            "specify both or neither of `prebuilt_geneseq` and `prebuilt_variants`"
        )
            
    pacbio_runs = pacbio_runs_from_config(config["pacbio_runs"])

    rule gene_sequence:
        """Get sequence of gene from PacBio amplicon."""
        input:
            gb=config["pacbio_amplicon"],
        output:
            codon=config["gene_sequence_codon"],
            prot=config["gene_sequence_protein"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "gene_sequence.txt"),
        script:
            "scripts/gene_sequence.py"


    rule align_parse_PacBio_ccs:
        """Align and parse PacBio CCS FASTQ file."""
        input:
            fastq=lambda wc: pacbio_runs.at[wc.pacbioRun, "fastq"],
            amplicon=config["pacbio_amplicon"],
            specs=config["pacbio_amplicon_specs"],
        output:
            outdir=directory(os.path.join(config["process_ccs_dir"], "{pacbioRun}")),
        conda:
            "environment_align_parse_PacBio_ccs.yml"
        log:
            os.path.join(config["logdir"], "align_parse_PacBio_ccs_{pacbioRun}.txt"),
        script:
            "scripts/align_parse_PacBio_ccs.py"


    rule analyze_pacbio_ccs:
        """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
        input:
            expand(rules.align_parse_PacBio_ccs.output.outdir, pacbioRun=pacbio_runs.index),
            config["pacbio_amplicon"],
            config["pacbio_amplicon_specs"],
            nb=os.path.join(config["pipeline_path"], "notebooks/analyze_pacbio_ccs.ipynb"),
        output:
            config["aligned_ccs_file"],
            nb="results/notebooks/analyze_pacbio_ccs.ipynb",
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "analyze_pacbio_ccs.txt"),
        shell:
            """
            papermill {input.nb} {output.nb} &> {log}
            head {output[0]}
            cat results/process_ccs/LibA_211105/readstats.csv
            """


    rule build_pacbio_consensus:
        """Build PacBio consensus sequences for barcodes."""
        input:
            config["aligned_ccs_file"],
            config["gene_sequence_codon"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/build_pacbio_consensus.ipynb"
            ),
        output:
            config["nt_variants"],
            nb="results/notebooks/build_pacbio_consensus.ipynb",
        params:
            config["max_ccs_error_rate"],
            config["consensus_params"],
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "build_pacbio_consensus.txt"),
        shell:
            "papermill {input.nb} {output.nb}"


    rule build_codon_variants:
        """Build codon-variant table."""
        input:
            config["nt_variants"],
            config["gene_sequence_codon"],
            config["gene_sequence_protein"],
            config["site_numbering_map"],
            config["mutation_design_classification"],
            config["neut_standard_barcodes"],
            nb=os.path.join(config["pipeline_path"], "notebooks/build_codon_variants.ipynb"),
        output:
            config["codon_variants"],
            nb="results/notebooks/build_codon_variants.ipynb",
        conda:
            "environment.yml"
        log:
            os.path.join(config["logdir"], "build_codon_variants.txt"),
        shell:
            "papermill {input.nb} {output.nb} &> {log}"
