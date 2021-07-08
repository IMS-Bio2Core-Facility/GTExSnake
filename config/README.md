# Configuration

> :warning: **Re-running a snakemake workflow**
> If you want to run this pipeline multpile times,
> then be sure to move or delete the final `results/process/sorted_isoforms.xlsx` file.
> Snakemake will not re-run a pipeline if it detects the output file.

There are currently five (5) points where the user can customise the pipeline:

## Number of Threads

Unless you are querying a huge number of genes,
I find 6 cores to be sufficient to keep things moving.
If you want to use more or less,
just change the value passed to `--cores`.
Each rule that implements concurrency will then use this many.
Remember, more isn't always faster.
On my laptop,
setting cores greater than the number of physical cores (not threads!)
in my machine gets me no improvement.
**YMMV**

## Genes of Interest

You will almost certainly have a different gene list than me!
If its short,
you can specify it as a parameter like so:

```shell
snakemake --use-conda --use-singularity --cores 6 --config gene_ids=\["GENE1","GENE2"\]
```

Alternatively, you can specify it as a list under the `gene_ids` parameter in the
configuration file located at `configuration/snakemake.yaml`.

## Region of Interest

The same goes for your region of interest.
It can be specified as the `region` parameter in `configuration/snakemake.yaml`,
or passed as below:

```shell
snakemake --use-conda --use-singularity --cores 6 --config region="Brain_Hypothalamus"
```

> :warning: Accepted Regions
> GTEx is picky about formatting,
> so be sure your region is from their approved list.
> You can find that list [here](HOLDING).

## Gencode URL

This is the URL to query for the Gencode reference.
We recommend **not** changing this unless you have a very strong reason.
Things start breaking rapidly if you use different refernces than GTEx.

## MANE URL

This is the URL to query for the MANE.
It *should* be up-to-date.
If we are late on an update,
or you would like to use either an older or new version,
change this.
