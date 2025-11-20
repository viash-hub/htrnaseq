# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

set -eo pipefail

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

nextflow \
  run . \
  -main-script src/workflows/runner/test.nf \
  -config ./src/config/labels.config \
  -entry test_wf \
  -resume \
  -profile docker,local

nextflow \
  run . \
  -main-script src/workflows/runner/test.nf \
  -config ./src/config/labels.config \
  -entry test_wf_with_lanes \
  -resume \
  -profile docker,local


nextflow \
  run . \
  -main-script src/workflows/runner/test.nf \
  -config ./src/config/labels.config \
  -entry test_wf_only_one_eset \
  -resume \
  -profile docker,local

nextflow \
  run . \
  -main-script src/workflows/runner/test.nf \
  -config ./src/config/labels.config \
  -entry test_same_experiment_id_different_projects_wf \
  -resume \
  -profile docker,local