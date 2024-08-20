# Viash-enabled HT-RNAseq pipeline

## Introduction
__TODO__: Add a description of the pipeline here.

## Test data

The pipeline can be run by creating a `params.yaml` file like this:

```yaml
param_list:
  - input_r1: "test_resources/100k/SRR14730301/VH02001612_S9_R1_001.fastq"
    input_r2: "test_resources/100k/SRR14730301/VH02001612_S9_R2_001.fastq"
    genomeDir: "test_resources/genomeDir/gencode.v41.star.sparse"
    barcodesFasta: "test_resources/2-wells.fasta"
    id: sample_one
  - input_r1: "test_resources/100k/SRR14730302/VH02001614_S8_R1_001.fastq"
    input_r2: "test_resources/100k/SRR14730302/VH02001614_S8_R2_001.fastq"
    genomeDir: "test_resources/genomeDir/gencode.v41.star.sparse"
    barcodesFasta: "test_resources/2-wells.fasta"
    id: sample_two
```

and then:

```bash
viash ns build --setup cb
nextflow run . -main-script target/nextflow/workflows/htrnaseq/main.nf \
  -profile docker \
  -c target/nextflow/workflows/htrnaseq/nextflow.config \
  -params-file params.yaml \
  -publish_dir output \
  -resume
```

Or, by running `src/workflows/htrnaseq/integration_test.sh`.
