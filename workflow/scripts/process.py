# -*- coding: utf-8 -*-
"""Process data to xlsx.

Here,
the input data from the previous three steps is combined and written to XLSX.
The GTEx data is merged BioMart data using an outer merge,
so as to keep all entries.
Then,
the MANE data is added using a left merge,
so as to only keep the data from the GTEx query.
"""
if __name__ == "__main__":
    import concurrent.futures

    import pandas as pd
    from gtexquery.data_handling.process import merge_data
    from gtexquery.logs.get_logger import get_logger

    INS = snakemake.input  # noqa: F821
    LOGS = snakemake.log[0]  # noqa: F821
    OUTS = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOGS)

    mane = pd.read_csv(INS["mane"], index_col=0)

    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
        ex.map(
            merge_data, INS["gtex"], INS["bm"], [mane] * len(INS["gtex"]), OUTS["data"]
        )
