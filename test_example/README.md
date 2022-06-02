# Test example of running pipeline
This subdirectory contains a test example of running the pipeline.
It is used to illustrate how the pipeline works and test it.

Note that this example differs from "normal" use as it is a subdirectory of the pipeline, whereas normally the pipeline is a subdirectory (submodule) of the analysis.
For this reason, the values of `pipeline_path` and `docs` in [config.yaml](config.yaml) should be changed in most other use cases as indicated in the comments in that file.
However, otherwise this test example can be used as a template.

The contents are:

 - [config.yaml](config.yaml): configuration for the `snakemake` pipeline.
 - [./results](results): results created by running the pipeline are here, only some are tracked in repo.
 - [Snakefile](Snakefile): top-level `snakemake` file.
 - [./data](data): input data
 - [./sequencing_data](sequencing_data): sequencing data for this test example; for your pipelines, the sequencing data probably will not be stored in repo.
 - [.gitignore](.gitignore): specifies which files to track, such as only tracking selected results.

You can run the pipeline with:

    snakemake --use-conda -j <n_jobs>

Note the pipeline uses the `conda` environment in [../environment.yml](../environment.yml), so if you have already built that environment (named `dms-vep-pipeline`), you can also just run:

    conda activate dms-vep-pipeline
    snakemake -j <n_jobs>

The sphinx HTML documentation is rendered to [https://dms-vep.github.io/dms-vep-pipeline/](https://dms-vep.github.io/dms-vep-pipeline/).
