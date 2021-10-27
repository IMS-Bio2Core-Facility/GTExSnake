rule process:
    input:
        gtex=rules.request.output.data,
        bm=rules.biomart.output.data,
        mane="resources/MANE.csv",
    output:
        data=report(
            expand("results/process/{gene}_isoforms.csv", gene=config["gene_ids"]),
            caption="../report/processed_data.rst",
        ),
    log:
        "results/logs/process.log",
    benchmark:
        "results/benchmarks/process.txt"
    threads: workflow.cores
    conda:
        "../envs/process.yaml"
    script:
        "../scripts/process.py"
