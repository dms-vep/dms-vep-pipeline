# Pipeline for analyzing deep mutational scanning (DMS) of viral entry proteins (VEPs)

## Details on setup of this repository

### Code formatting and linting
The Python code should be formatted with [black](https://black.readthedocs.io/) by running `black .`

Comparable formatting is done for the `Snakefile` with [snakefmt](https://github.com/snakemake/snakefmt) by running `snakefmt .`.

The `Snakefile` should be linted with `cd test_example; snakemake --lint; cd ..` (note that this lints starting with the top `Snakefile` that runs the entire test pipeline).

The code and Jupyter notebooks are linted with [flake8_nb](https://flake8-nb.readthedocs.io/) by running `flake8_nb`.

### Testing of pipeline on Travis CI
The pipeline is tested on Travis CI by checking all the formatting and linting above, and then also runing the pipeline on the example in [./test_example](test_example).

### Stripping of Jupyter notebook output
The repo was configured to strip output from Jupyter notebooks as described [here](http://mateos.io/blog/jupyter-notebook-in-git/) by running:

    nbstripout --install --attributes .gitattributes

### Tracking of large files with `git lfs`
The large data files for the test example in [./test_example/sequencing_data/](test_example/sequencing_data/) are tracked with [git lfs](http://arfc.github.io/manual/guides/git-lfs).
Note the first time you set up the repo, you have to run:

    git lfs install
