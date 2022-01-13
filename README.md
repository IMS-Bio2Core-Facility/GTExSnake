# Snakemake Workflow: GTExSnake

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CI/CD](https://github.com/IMS-Bio2Core-Facility/GTExSnake/actions/workflows/cicd.yaml/badge.svg)](https://github.com/IMS-Bio2Core-Facility/GTExSnake/actions/workflows/cicd.yaml)
[![Codestyle: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Codestyle: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

A fully concurrent pipeline for querying transcript-level GTEx data in specific tissues

**_Find us in the "Standardized Usage" Section of the [Snakemake Workflow Catalog][sm_wc]_**

If you find the project useful,
leaves us a star on [github][stars]!

A list of changes can be found in our [CHANGELOG](./CHANGELOG.md).

## Motivation

There are a number of circumstances where transcript level expressed data for a
specific tissue is highly valuable.
For tissue-dependent expression data,
there are few resources better than GTEx.
In this case, the `medianTranscriptExpression` query provides the necessary data.
It returns the median expression of each transcript for a gene in a given tissue.
Here, we query a list of genes against a region-specific subset of GTEx.

## Notes on Installation

A full walkthrough on how to install and use this pipeline
can be found [here][sm_wc].

To take advantage of Singularity,
you'll need to install that separately.
If you are running on a Linux system,
then singularity can be installed from conda like so:

```shell
conda install -n snakemake -c conda-forge singularity
```

It's a bit more challenging for other operating systems.
Your best bet is to follow their instructions
[here][sing_install].
But don't worry!
**Singularity is _not_ regquired!**
Snakemake will still run each step in its own Conda environment,
it just won't put each Conda environment in a container.

### Get the Source Code

Alternatively,
you may grab the source code.
You likely will not need these steps if you aren't planning to contribute.

Navigate to our [release][releases]
page on github and download the most recent version.
The following will do the trick:

```shell
curl -s https://api.github.com/repos/IMS-Bio2Core-Facility/GTExSnake/releases/latest |
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
git clone https://github.com/IMS-Bio2Core-Facility/GTExSnake
```

> :warning: **Heads Up!**
> The bleeding edge may not be stable,
> as it contains all active development.

However you choose to install it,
`cd` into the directory.

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
it will be available on the [Snakemake workflow catalog][sm_wc]
in addition to the usual availablility on [Github][repo].

Unfortunately,
any query to an unstable API is inherently not reproducible.
Thus,
changes in BioMart or GTEx could impact the results.
We recognise this as an inherent limitation,
and will do our best to keep abreast of API changes that impact the pipeline.

## Data

The pipeline requires no input data other than a list of gene names specified in
[config/config.yaml](./config/config.yaml).

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

## Contributing

If you are interested in helping us improve the pipeline,
pleare see our guides on [contributing](./CONTRIBUTING.md)
and be sure to abide by our [code of conduct](./CODE_OF_CONDUCT.md)!

[stars]: https://github.com/IMS-Bio2Core-Facility/GTExSnake/stargazers "Stargazers"
[conda]: https://docs.conda.io/en/latest/ "Conda"
[snakemake]: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html "Snakemake"
[singularity]: https://sylabs.io/singularity/ "Singularity"
[mambaforge]: https://github.com/conda-forge/miniforge#mambaforge "Mambaforge"
[sing_install]: https://sylabs.io/guides/3.8/admin-guide/installation.html#installation-on-windows-or-mac "Singularity Install"
[releases]: https://github.com/IMS-Bio2Core-Facility/GTExSnake/releases "GTExSnake Releases"
[sm_wc]: https://github.com/IMS-Bio2Core-Facility/GTExSnake "Usage Instructions"
[repo]: https://github.com/IMS-Bio2Core-Facility/GTExSnake "GTExSnake Repository"
[MANE]: https://www.ncbi.nlm.nih.gov/refseq/MANE/ "MANE"
