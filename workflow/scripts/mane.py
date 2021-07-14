# -*- coding: utf-8 -*-
"""Download the summary file from the MANE project.

Matched annotation from NCBI and EMBL-EBI, or `MANE`_,
is an initiative to "provide a minimal set of matching RefSeq and Ensembl transcripts
of human protien-coding genes."
Though it represents only a minimal set and has not yet achieved full genome coverage,
it is an incredibly valuable resource for matching these two references together.

Here, the ``MANE.GRCh38.v0.93.summary.txt.gz`` file is retrieved with the ever useful
``pandas.read_csv``. After renaming several columns for compatibility with GTEx and
BioMart data, the verion numbers are stripped from the ENSG, ENST, and RefSeq IDs,
as GTEx does not use the most recent Gencode.

.. _MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/
"""
if __name__ == "__main__":

    import pandas as pd

    mane = pd.read_csv(
        "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.93/MANE.GRCh38.v0.93.summary.txt.gz",
        sep="\t",
        header=0,
    )
    mane = mane.rename(
        columns={
            "Ensembl_Gene": "gencodeId",
            "symbol": "geneSymbol",
            "Ensembl_nuc": "transcriptId",
            "RefSeq_nuc": "refseq",
        }
    )

    mane.loc[:, ["gencodeId", "transcriptId", "refseq"]] = mane.loc[
        :, ["gencodeId", "transcriptId", "refseq"]
    ].apply(lambda x: x.str.split(".").str.get(0))

    mane.to_csv("resources/MANE.csv")
