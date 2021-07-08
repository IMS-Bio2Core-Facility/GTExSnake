# Snakemake Workflow: GTExSnake

> :warning: **This workflow will not currently run**
> The home repository [here][home]
> is currently being refactored to split the project into 2 repositories:
>
> 1) This one, containing the pipeline.
> 2) A forthcoming one, to be a pip package containing the code.

A fully concurrent pipeline for querying transcript-level GTEx data in specific tissues

If you find the project useful,
leaves us a star on [github][stars]!

[home]: https://github.com/IMS-Bio2Core-Facility/BIC086 "Original repository"
[stars]: HOLDING "Stargazers"

## Motivation

There are a number of circumstances where transcript level expressed data for a
specific tissue is highly valuable.
For tissue-dependent expression data,
there are few resources better than GTEx.
In this case, the `medianTranscriptExpression` query provides the necessary data.
It returns the median expression of each transcript for a gene in a given tissue.
Here, we query a list of genes against a region-specific subset of GTEx.

## Pipeline

### Necessary Software

This pipeline needs [conda][conda]
and [snakemake][snakemake]
installed,
and runs best if you also have [singularity][singularity]
installed,
though it's not required.

Snakemake recommends using [mambaforge][mambaforge]
as your base conda,
which I would also recommend.
Installation instructions are at the above link.
If you prefer a vanilla conda installation,
you can always try `mamba` following the instructions at the above snakemake link.
Once you have conda installed,
install snakemake as outlined on their page
(again, see the above link).

Singularity can be a bit of a nuisance to install
(at least in my experience).
I would recommend a Linux OS
(it's most straightforward on CentOS/RHEL  - hello AWS!)
and following the instructions in their repository
[here][sing_install].
If you decide not to install Singularity,
don't worry!
Snakemake will still run each step in its own Conda environment,
it just won't put each Conda environment in a container.

[conda]: https://docs.conda.io/en/latest/ "Conda"
[snakemake]: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html "Snakemake"
[singularity]: https://sylabs.io/singularity/ "Singularity"
[mambaforge]: https://github.com/conda-forge/miniforge#mambaforge "Mambaforge"
[sing_install]: https://github.com/sylabs/singularity/blob/master/INSTALL.md "Singularity Install"

### Get the Source Code

Navigate to our [release](HOLDING)
page on github and download the most recent version.
The following will do the trick:

```shell
curl -s https://api.github.com/repos/IMS-Bio2Core-Facility/BIC086/releases/latest |
grep tarball_url |
cut -d " " -f 4 |
tr -d '",' |
xargs -n1 curl -sL |
tar xzf -
```

After querying the github api to get the most recent release information,
we grep for the desired URL,
split the line and extract the field,
trim superfluous characters,
use `xargs` to pipe this to `curl` while allowing for re-directs,
and un-tar the files.
Easy!

Alternatively,
for the bleeding edge,
please clone the repo like so:

```shell
git clone https://github.com/IMS-Bio2Core-Facility/BIC086
```

> :warning: **Heads Up!**
> The bleeding edge may not be stable,
> as it contains all active development.

However you choose to install it,
`cd` into the directory and activate your snakemake conda environment.

### Running

Once you've installed the above software,
and fetched the code,
running the pipeline is as simple as:

```shell
snakemake --use-conda --use-singularity --cores 6
```

If you aren't using `singularity`,
then leave off the apropriate flag, as so:

```shell
snakemake --use-conda --cores 6
```

And `snakemake` will automatically leave it off.

### Configuration

For more on how to configure your custom run,
please see [config/README.md][config].

[config]: config/README.md

## Reproducibility

Reproducibility results are a cornerstone of the scientific process.
By running the pipeline with `snakemake` in a `docker` image using `conda` environments,
we ensure that no aspect of the pipeline is left to chance.
You will get our analysis,
as we ran it,
with the software versions,
as we used them.

We also strive to make this pipeline as FAIR/O compliant as possible.
To that end,
it is available both on the [Snakemake workflow catalog][snake_catalog]
and [Workflowhub][workflowhub],
in addition to the usual availablility on [Github][repo].

Unfortunately,
any query to an unstable API is inherently not reproducible.
Thus,
changes in BioMart or GTEx could impact the results.
We recognise this as an inherent limitation,
and will do our best to keep abreast of API changes that impact the pipeline.

[snake_catalog]: HOLDING
[workflowhub]: HOLDING
[repo]: HOLDING

## Data

The pipeline requires no input data other than a list of gene names specified in
`configuration/snakemake.yaml`.

## References

It is surprisingly challenging to align RefSeq IDs and Ensembl IDs.
This is further complicated because GTEx uses Gencode26 under the hood.
As this is not the most up-to-date version,
it actually proved quite frustrating to find the desired version numbers for each gene.
To combat this,
this pipeline takes 3 different approaches in parallel:

1. Gencode v26 GTF annotations are downloaded from EBI,
   so the user only needs to supply gene names.
1. A query is made to BioMart to retrieve RefSeq IDs for each ENST returned by GTEx.
1. Data from [MANE][mane] is added to help identify consensus transcripts.

[MANE]: https://www.ncbi.nlm.nih.gov/refseq/MANE/ "MANE"

## Contributing

If you are interested in helping us improve the pipeline,
pleare see our guides on [contributing][contributing]
and be sure to abide by our [code of conduct][conduct]!

[contributing]: CONTRIBUTING.md
