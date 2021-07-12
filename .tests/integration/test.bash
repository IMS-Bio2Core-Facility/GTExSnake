#!/bin/usr/env bash

# All data is either retrieved by a step in the pipeline
# or it is contained in the configuration
# so the only thing necessary for the integration tests is...
snakemake -j6 --use-conda --use-singularity
