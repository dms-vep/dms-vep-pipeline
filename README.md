# Pipeline for analyzing deep mutational scanning (DMS) of viral entry proteins (VEPs)
[![Build Status](https://github.com/dms-vep/dms-vep-pipeline/actions/workflows/test.yaml/badge.svg)](https://github.com/dms-vep/dms-vep-pipeline/actions/workflows/test.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


## Overview
This repository contains a [snakemake](https://snakemake.readthedocs.io/) pipeline for analysis of deep mutational scanning of barcoded viral entry proteins.
In order to use the pipeline, you can include this repository as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in your own repo, which will contain the data and project-specific code.

In other words, if the repository for your specific deep mutational scanning project is called `<my_dms_repo>`, you would add [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) as a submodule to that, while the master [snakemake](https://snakemake.readthedocs.io/) `Snakefile`, its configuration (`config.yaml`), its input data, etc would reside in `<my_dms_repo>`.
In other words, the directory structure would look like this:

```
<my_dms_repo>
├── dms-vep-pipeline [added as git submodule]
├── README.md [README for main project]
├── Snakefile [top-level snakemake file]
├── config.yaml [configuration for Snakefile]
├── data [subdirectory with input data]
├── results [subdirectory with results created by pipeline]
├── docs [sphinx summary of results created by pipeline]
└── <other files / subdirectories that are part of project>
```

The top-level `Snakefile` then [includes](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes) the [snakemake](https://snakemake.readthedocs.io/) rules defined by the `*.smk` files in the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) submodule, and uses them to run the analysis.
This also requires properly setting up the top-level `config.yaml` to specify details for your project (see more below).

## Test example
A test example use of the pipeline is in [./test_example](test_example).
In order to make the example contained in the pipeline, the organization is a bit different for this [./test_example](test_example): it is contained as a subdirectory of the pipeline whereas for actual use of the repo you will make the pipeline a submodule of `<my_dms_repo>` as described above.
Therefore, your `config.yaml` will have different values for `pipeline_path` and `docs` as indicated in the comments in [./test_example/config.yaml](test_example/config.yaml).

Despite these differences, [./test_example](test_example) provides an example of how to set up your repo.
Running the `Snakefile` in [./test_example](test_example) performs the whole analysis for the test example and creates the sphinx rendering in [./docs](docs) which can be displayed via GitHub pages as here: [https://dms-vep.github.io/dms-vep-pipeline/](https://dms-vep.github.io/dms-vep-pipeline/).

## Running the pipeline
The `Snakefile` you create will include [pipeline.smk](pipeline.smk) (which has the analysis pipeline) and [docs.smk](docs.smk) (which builds the [sphinx](https://www.sphinx-doc.org/) HTML documentation).
You can also optionally add other rules into your `Snakefile`.
If they define an `output` named `nb` that is a Jupyter notebook (like some of the rules in [pipeline.smk](pipeline.smk) and its included `.smk` files), then that will be included into the HTML documentation.

You then run the pipeline with:

    snakemake -j <n_jobs> --use-conda

Or if you are only using the `dms-vep-pipeline` [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml) and have already built that, you can also just do:

    conda activate dms-vep-pipeline
    snakemake -j <n_jobs>

If the `./docs` output directory has already been built and you want to force a re-run, just delete it and then run above.

This will create the results in `./results/` and the HTML documentation in `./docs/`.
To display the HTML documentation via GitHub pages, set up your repo to serve documentation via GitHub pages from the `/docs` folder of the main (or master) branch [as described here](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).
The documentation will then be at `https://dms-vep.github.io/<my_dms_repo>` (assuming you are using the [https://github.com/dms-vep](https://github.com/dms-vep) organization; otherwise replace `dms-vep` with whatever account contains your repo).

Note that [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) has its own [conda](https://docs.conda.io/) environment specified in [environment.yml](environment.yml).
There is a separate environment, [environment_align_parse_PacBio_ccs.yml](environment_align_parse_PacBio_ccs.yml), for aligning and parsing the PacBio CCSs so that isn't re-run every time main environment is updated.

## Setting up the pipeline as a submodule
To add [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in your repo (`<my_dms_repo>`), do as follows.

        git submodule add https://github.com/dms-vep/dms-vep-pipeline

This adds the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline](dms-vep-pipeline), which can then be committed with:

    git commit -m 'added `dms-vep-pipeline` as submodule'

Note that if you want a specific commit or tag of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), follow the [steps here](https://stackoverflow.com/a/10916398):

    cd dms-vep-pipeline
    git checkout <commit>

and then `cd ../` back to the top-level directory of `<my_dms_repo>`, and add and commit the updated `dms-vep-pipeline` submodule.

You can also make changes to the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) submodule in `<my_dms_repo>` by going into that directory, making changes on a branch, and then pushing back to [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) and opening a pull request.

## Details on setup of this pipeline repository

### Structure of this repo
Here are the different contents of this repo:

 - [conf.py](conf.py): [sphinx](https://www.sphinx-doc.org/) configuration used to build HTML docs.
 - [docs.smk](docs.smk): [snakemake](https://snakemake.readthedocs.io/) rules used to build HTML docs.
 - [funcs.smk](funcs.smk) functions used in [snakemake](https://snakemake.readthedocs.io/) pipeline.
 - [LICENSE.txt](LICENSE.txt): license for pipeline.
 - [pipeline.smk](pipeline.smk): [snakemake](https://snakemake.readthedocs.io/) rules that run pipeline.
 - [build_variants.smk](build_variants.smk): [snakemake](https://snakemake.readthedocs.io/) rules for pipeline specific to building variants, included in [pipeline.smk](pipeline.smk).
 - [./scripts](scripts): subdirectory with scripts used by [snakemake](https://snakemake.readthedocs.io/) rules.
 - [./docs](docs): HTML documentation for the test example.
 - [environment.yml](environment.yml): the [conda](https://docs.conda.io/) environment for the pipeline.
 - [./notebooks](notebooks): Jupyter notebooks used by [snakemake](https://snakemake.readthedocs.io/) rules.
 - [README.md](README.md): README file describing pipeline.
 - [./test_example](test_example): test example that runs pipeline.
 - [.flake8_nb](.flake8_nb): configuration for [flake8_nb](https://flake8-nb.readthedocs.io/) linting.
 - [.gitattributes](.gitattributes)
 - [.gitignore](.gitignore): files to ignore in GitHub tracking
 - [.github/workflows/test.yaml](.github/workflows/test.yaml): specifying how testing on GitHub Actions is done.

### `conda` environment
The [conda](https://docs.conda.io/) environment for the pipeline is in [environment.yml](environment.yml).

### Code formatting and linting
The Python code should be formatted with [black](https://black.readthedocs.io/) by running `black .`

Comparable formatting is done for the `snakemake` file (`*.smk` files) with [snakefmt](https://github.com/snakemake/snakefmt) by running `snakefmt .`.
The overall `snakemake` pipeline is linted by going to [./test_example](test_example) and running `snakemake --lint`.

The code and Jupyter notebooks are linted with [flake8_nb](https://flake8-nb.readthedocs.io/) by running `flake8_nb`.

### Testing of pipeline with GitHub Actions
The pipeline is tested with GitHub Actions by checking all the formatting and linting above, and then also running the pipeline on the example in [./test_example](test_example).
See [.github/workflows/test.yaml](.github/workflows/test.yaml) for details.

### Stripping of Jupyter notebook output
The repo was configured to strip output from Jupyter notebooks as described [here](http://mateos.io/blog/jupyter-notebook-in-git/) by running:

    nbstripout --install --attributes .gitattributes

### Tracking of large files with `git lfs`
The large data files for the test example in [./test_example/sequencing_data/](test_example/sequencing_data/) are tracked with [git lfs](http://arfc.github.io/manual/guides/git-lfs).
Note the first time you set up the repo, you have to run:

    git lfs install

### Notes
If you add noteboooks in subindices, to avoid ``sphinx`` errors they must have the "orphan" tag: https://nbsphinx.readthedocs.io/en/0.8.8/orphan.html
