"""Create sub-index ``*.rst`` file for sphinx docs."""

import os
import sys


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

title = snakemake.params.title
subindex = os.path.splitext(os.path.basename(snakemake.output.subindex))[0]
to_replace = f"{subindex}_"
nblinks = [os.path.splitext(os.path.basename(nblink))[0] for nblink in snakemake.input]
analysis_nbs = "\n".join(
    f"- :doc:`{nblink.replace(to_replace, '')} <{nblink}>`" for nblink in nblinks
)

with open(snakemake.output.subindex, "w") as f_obj:
    f_obj.write(
        f"""\

{title}
{"-" * len(title)}

{analysis_nbs}

"""
    )
