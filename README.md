# Pipeline for analyzing deep mutational scanning (DMS) of viral entry proteins (VEPs)

## Details on setup of this repository

### Code formatting and linting
The code should be formatted with [black](https://black.readthedocs.io/) by running `black .`

### Stripping of Jupyter notebook output
The repo was configured to strip output from Jupyter notebooks as described [here](http://mateos.io/blog/jupyter-notebook-in-git/) by running:

    nbstripout --install --attributes .gitattributes
