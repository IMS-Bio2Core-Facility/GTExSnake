rule request:
    input:
        lut="resources/gencode.v26.annotation",
    params:
        gene_ids=config["gene_ids"],
        region=config["region"],
    output:
        data=expand("results/request/{gene}_message.csv", gene=config["gene_ids"]),
    log:
        "results/logs/request.log",
    benchmark:
        "results/benchmarks/request.txt"
    threads: workflow.cores
    conda:
        "../envs/request.yaml"
    script:
        "../scripts/request.py"


rule biomart:
    input:
        data=rules.request.output.data,
    output:
        data=expand("results/biomart/{gene}_message.csv", gene=config["gene_ids"]),
    log:
        "results/logs/biomart.log",
    benchmark:
        "results/benchmarks/biomart.txt"
    threads: workflow.cores
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/biomart.py"
