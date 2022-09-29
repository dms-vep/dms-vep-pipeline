"""Create ``index.rst`` file for sphinx docs."""

import os
import sys


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

github_user = snakemake.config["github_user"]
github_repo = snakemake.config["github_repo"]
github_branch = snakemake.config["github_branch"]
authors = snakemake.config["authors"]

github_url = f"https://github.com/{github_user}/{github_repo}"
blob_path = f"{github_url}/blob/{github_branch}"

docs_source_relpath = snakemake.params.docs_source_relpath
rulegraph = os.path.join(docs_source_relpath, snakemake.input.rulegraph)
filegraph = os.path.join(docs_source_relpath, snakemake.input.filegraph)

nbs_for_index = snakemake.params.nbs_for_index
analysis_nbs = "\n   ".join(nbs_for_index)

results_relpath = snakemake.params.results_relpath
data_file_links = []
for label, link in snakemake.params.data_files.items():
    if isinstance(link, tuple):
        assert len(link) == 2
        link = link[1]
    assert isinstance(link, str)
    data_file_links.append(f"- `{label} <{blob_path}/{results_relpath}/{link}>`_")
data_file_links = "\n".join(data_file_links)

with open(snakemake.output.index, "w") as f_obj:
    f_obj.write(
        f"""\
{snakemake.config["description"]}
{"=" * len(snakemake.config["description"])}

This page documents the data analysis.
For the actual code, see {github_url}

Study by {authors}.

Workflow
--------
Below is the rulegraph for the `snakemake <https://snakemake.readthedocs.io/>`_ workflow.
Click :download:`here <{filegraph}>` for the more detailed filegraph.

.. image:: {rulegraph}
   :width: 800

Analysis notebooks
------------------
Many of the plots in these notebooks are interactive, so try mousing over points for
details, using dropdowns, etc.

.. toctree::
   :maxdepth: 1

   {analysis_nbs}

Data files
----------
{data_file_links}

"""
    )

    # link HTML plots: https://stackoverflow.com/a/67997311
    if snakemake.params.have_func_selections:
        f_obj.write(
            f"""\
Interactive plots of mutation functional effects
------------------------------------------------
- `Observed phenotype effects <{os.path.basename(snakemake.input.muteffects_observed)}>`_
- `Latent phenotype effects <{os.path.basename(snakemake.input.muteffects_latent)}>`_

"""
        )

    if snakemake.params.have_antibody_selections:
        f_obj.write(
            """\
Interactive plots of mutation antibody escape
---------------------------------------------
"""
        )
        for html_plot in sorted(snakemake.input.antibody_escape_plots):
            base_plot = os.path.basename(html_plot)
            label = (
                os.path.splitext(base_plot)[0]
                .replace("_formatted", " ")
                .replace("_", " ")
            )
            f_obj.write(f"- `{label} <{base_plot}>`_\n")
