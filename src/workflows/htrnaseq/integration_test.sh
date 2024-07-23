#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# Make sure the workflow is built
viash ns build -q workflows

export NXF_VER=24.04.3

nextflow run . \
  -main-script target/nextflow/workflows/htrnaseq/main.nf \
  -params-file params1.yaml \
  -with-docker \
  --publish_dir output \
  -resume

# bin/nextflow run . \
#   -main-script src/workflows/htrnaseq_wf/test.nf \
#   -entry test_wf_NextSeq550 \
#   -profile docker,local \
#   -resume \
#   --mappingReferenceRoot testData/genomeDir/subset \
#   --barcodesReferenceRoot testData \
#   --compoundAnnotationRoot testData \
#   -ansi-log false \
#   --publish_dir htrnaseq_test_results
#
