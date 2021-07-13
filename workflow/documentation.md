# Technical Documentation

This information is also available as the docstring for the analysis script,
where apropriate.

## Rule: ids

Retrieves the Gencode v26 annotations.

The website from which the annotations are retrieved is specified in
[config/config.yaml](../config/config.yaml)
under `gencode_url` and can be altered by the user.

This rule retrieves the annotated gtf file using `wget`
and decompresses it using `gzip`.
The apropriate fields are selected using `cut`
and unnecessary characters (here, ; and ") are removed with `tr`
before the file is directed to `resources/gencode.v26.annotation`.

## Rule: mane

Download the summary file from the MANE project.

The website from which the annotations are retrieved is specified in
[config/config.yaml](../config/config.yaml)
under `mane_url` and can be altered by the user.

These comments are taken from the `workflow/scripts/mane.py` docstring:

> Matched annotation from NCBI and EMBL-EBI, or [MANE][mane],
> is an initiative to "provide a minimal set of matching RefSeq and Ensembl transcripts
> of human protien-coding genes."
> Though it represents only a minimal set and has not yet achieved full genome coverage,
> it is an incredibly valuable resource for matching these two references together.
>
> Here, the `MANE.GRCh38.v0.93.summary.txt.gz` file is retrieved with the ever useful
> `pandas.read_csv`. After renaming several columns for compatibility with GTEx and
> BioMart data, the verion numbers are stripped from the ENSG, ENST, and RefSeq IDs,
> as GTEx does not use the most recent Gencode.

## Rule: request

Make a GET request to GTEx medianTranscriptExpression.

The brain region and gene ids for which the request is made are specified in
[config/config.yaml][config]
under `region` and `gene_ids`,
respectively,
and can be altered by the user.

Concurrency is achieved by using a `concurrent.futures.ThreadPoolExecutor`
to map the `gtex_reguest` function from the
[GTExQuery][gtexquery]
package over the region, genes, and output files.

These comments are taken from `workflow/scripts/request.py` docstring:

> This step queries the GTEx API for transcript expression data in the region
> specified by the user,
> using a user provided list of genes names.
> To simplify the use experience,
> the user should provide common gene names.
> These are then automatically converted to the GTEx-required Ensembl IDs
> by referencing Gencode v26 as this is the version used by GTEx.
> If a gene name is not found in Gencode,
> then a warning is dumped to the logs,
> and a blank file created to propagate this error downstream.
>
> As the GTEx API is quite straightforward,
> these queries can be made using the standard [requests.session][session] object.
> Data were pulled from the `gtex_v8` dataset limited to the
> region specified by the user.
>
> There are a number of bugs in the GTEx API that necesitated some work arounds.
> The data are returned in a TSV format,
> as attempts to return in a CSV format produced a tab-delimited first line
> and comma-delimited body.
> Additionally,
> the `hclust` option cannot be passed as it is force-set to `true`,
> no matter the user-passed value.
> This is circumvented by not passing the value in the first place.

## Rule: biomart

Use BioMart to convert ENST IDs to RefSeq IDs.

There are no user specified parameters for this step.

Concurrency is achieved by using a `concurrent.futures.ThreadPoolExecutor`
to map the `concurrent_biomart` function from the
[GTExQuery][gtexquery]
package over the input and output files.

These comments are taken from `workflow/scripts/biomart.py` docstring:

> Once the data has been retrieved from GTEx,
> the RefSeq IDs need to be added.
> The most straightforward way to achieve this is to query BioMart,
> using the ENST as a filter.
>
> This is achieved using the `BioMart` class from the [bioservices][bioservices] package.
> The query is made against the `hsapiens_gene_ensembl` dataset,
> filtered by ENST,
> and the ENSG ID, ENST ID, HGNC symbol, and RefSeq ID are returned.
>
> Due to the differing standards used by Ensembl and RefSeq,
> not all ENST IDs have a corresponding RefSeq ID and,
> confusingly,
> some ENST IDs map to multiple RefSeq IDs.
> It was the latter case that necesitated the addition of data from MANE
> to provide some clarity to the situation.

## Rule: process

Process data to XLSX.

There are no user specified parameters for this step.

Concurrency is achieved using the `Pipeline` class from the
[GTExQuery][gtexquery]
package to map the file reading, merging, and writing over the
input and output files.

These comments are taken from `workflow/scripts/process.py` docstring:

> Here,
> the input data from the previous three steps is combined and written to XLSX.
> The GTEx data is merged BioMart data using an outer merge,
> so as to keep all entries.
> Then,
> the MANE data is added using a left merge,
> so as to only keep the data from the GTEx query.

[mane]: https://www.ncbi.nlm.nih.gov/refseq/MANE/ "MANE"
[session]: https://docs.python-requests.org/en/master/user/advanced/ "Requests Package"
[bioservices]: https://github.com/cokelaer/bioservices "Bioservices Package"
[gtexquery]: https://github.com/IMS-Bio2Core-Facility/GTExQuery "GTExQuery"
