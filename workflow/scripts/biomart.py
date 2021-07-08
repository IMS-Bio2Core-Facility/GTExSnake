# -*- coding: utf-8 -*-
"""Use BioMart to convert ENST IDs to RefSeq IDs.

Once the data has been retrieved from GTEx,
the RefSeq IDs need to be added.
The most straightforward way to achieve this is to query BioMart,
using the ENST as a filter.

This is achieved using the ``BioMart`` class from the `bioservices`_ package.
The query is made against the ``hsapiens_gene_ensembl`` dataset,
filtered by ENST,
and the ENSG ID, ENST ID, HGNC symbol, and RefSeq ID are returned.

Due to the differing standards used by Ensembl and RefSeq,
not all ENST IDs have a corresponding RefSeq ID and,
confusingly,
some ENST IDs map to multiple RefSeq IDs.
It was the latter case that necesitated the addition of data from MANE
to provide some clarity to the situation.

.. _bioservices: https://github.com/cokelaer/bioservices
"""
if __name__ == "__main__":
    import concurrent.futures

    # import pandas as pd
    from logs.get_logger import get_logger
    from multithreading.biomart import concurrent_biomart

    INPUTS = snakemake.input  # noqa: F821
    LOGS = snakemake.log[0]  # noqa: F821
    OUTS = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOGS)

    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
        ex.map(concurrent_biomart, INPUTS["data"], OUTS["data"])
