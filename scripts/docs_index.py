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
dag = os.path.join(docs_source_relpath, snakemake.input.dag)

nbs_for_index = snakemake.params.nbs_for_index
analysis_nbs = "\n   ".join(nbs_for_index)

results_relpath = snakemake.params.results_relpath
data_file_links = "\n".join(
    [
        f"- `{label} <{blob_path}/{results_relpath}/{link}>`_"
        for label, link in snakemake.params.data_files.items()
    ]
)

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
Click :download:`here <{filegraph}>` for the more detailed filegraph,
and :download:`here <{dag}>` for the even more detailed DAG (directed acyclic graph).

.. image:: {rulegraph}
   :width: 600

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
