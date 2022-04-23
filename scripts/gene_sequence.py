"""Implements ``snakemake`` rule to get gene sequence."""


import sys

import Bio.SeqIO


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

print(f"Reading amplicon from {snakemake.input.gb=}")
amplicon = Bio.SeqIO.read(snakemake.input.gb, "genbank")

gene_feature = [f for f in amplicon.features if f.type == "gene"]
if len(gene_feature) != 1:
    raise ValueError("failed to find exactly one feature of type 'gene'")
geneseq = gene_feature[0].extract(amplicon).seq
print(f"Found gene of length {len(geneseq)}")

# first make sure gene can translate as CDS
_ = geneseq.translate(cds=True)
# now translate without cds flag to keep any stop codons
protseq = geneseq.translate()

for outfile, seq in [
    (snakemake.output.codon, geneseq),
    (snakemake.output.prot, protseq),
]:
    print(f"Writing to {outfile}")
    with open(outfile, "w") as f:
        f.write(f">gene\n{str(seq)}\n")
