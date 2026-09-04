#!/usr/bin/env bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q well_demultiplex_runner --setup cb

nextflow run . \
  -main-script src/workflows/well_demultiplex_runner/test.nf \
  -entry test_wf \
  -profile docker,local \
  -c src/config/labels.config \
  -resume