#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# Make sure the workflow is built
viash ns build --setup cb --parallel

export NXF_VER=24.04.4

nextflow \
  run . \
  -main-script src/workflows/well_demultiplex/test.nf \
  -config ./src/config/labels.config \
  -entry test_wf \
  -resume \
  -profile docker,local \
  --publish_dir output \
