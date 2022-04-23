# Pipeline for analyzing deep mutational scanning (DMS) of viral entry proteins (VEPs)

## Details on setup of this repository

### Code formatting and linting
The Python code should be formatted with [black](https://black.readthedocs.io/) by running `black .`

Comparable formatting is done for the `Snakefile` with [snakefmt](https://github.com/snakemake/snakefmt) by running `snakefmt .`.

The `Snakefile` should be linted with `cd test_example; snakemake --lint; cd ..` (note that this lints starting with the top `Snakefile` that runs the entire test pipeline).

The code and Jupyter notebooks are linted with [flake8_nb](https://flake8-nb.readthedocs.io/) by running `flake8_nb`.

### Stripping of Jupyter notebook output
The repo was configured to strip output from Jupyter notebooks as described [here](http://mateos.io/blog/jupyter-notebook-in-git/) by running:

    nbstripout --install --attributes .gitattributes
