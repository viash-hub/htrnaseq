# HT-RNAseq - A pipeline for processing high-throughput RNA-seq data

## Introduction
__TODO__: Add a description of the pipeline here.

## Test data

As test data, we use [a DRUGseq dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176150) from the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).

The original data has been (partly) subsampled to reduce the test runtime. We used [seqtk](https://github.com/lh3/seqtk) for this with a seed of 1, e.g.:

```bash
seqtk sample -s1 orig/SRR14730302/VH02001614_S8_R1_001.fastq.gz 10000 > 10k/SRR14730302/VH02001614_S8_R1_001.fastq.gz
```

The data is available at: `gs://viash-hub-test-data/htrnaseq/v1/`:

```
❯ gcstree -f viash-hub-test-data/htrnaseq/v1/
viash-hub-test-data
└── htrnaseq
    └── v1
        ├── [  48]  2-wells.fasta
        ├── [465.3K]  GSE176150_metadata.csv
        ├── 100k
        │   ├── SRR14730301
        │   │   ├── [8.5M]  VH02001612_S9_R1_001.fastq
        │   │   └── [14.9M]  VH02001612_S9_R2_001.fastq
        │   └── SRR14730302
        │       ├── [8.5M]  VH02001614_S8_R1_001.fastq.gz
        │       └── [14.9M]  VH02001614_S8_R2_001.fastq.gz
        ├── 10k
        │   ├── SRR14730301
        │   │   ├── [845.4K]  VH02001612_S9_R1_001.fastq
        │   │   └── [1.5M]  VH02001612_S9_R2_001.fastq
        │   └── SRR14730302
        │       ├── [845.3K]  VH02001614_S8_R1_001.fastq.gz
        │       └── [1.5M]  VH02001614_S8_R2_001.fastq.gz
        └── orig
            ├── [20.4G]  SRR14730301
            │   └── [20.4G]  SRR14730301
            ├── SRR14730301
            │   ├── [9.1G]  VH02001612_S9_R1_001.fastq.gz
            │   └── [22.0G]  VH02001612_S9_R2_001.fastq.gz
            ├── [16.9G]  SRR14730302
            │   └── [16.9G]  SRR14730302
            ├── SRR14730302
            │   ├── [7.6G]  VH02001614_S8_R1_001.fastq.gz
            │   └── [18.0G]  VH02001614_S8_R2_001.fastq.gz
            ├── [18.0G]  SRR14730303
            │   └── [18.0G]  SRR14730303
            ├── SRR14730303
            │   ├── [8.1G]  VH02001618_S7_R1_001.fastq.gz
            │   └── [19.2G]  VH02001618_S7_R2_001.fastq.gz
            ├── [16.5G]  SRR14730304
            │   └── [16.5G]  SRR14730304
            ├── SRR14730304
            │   ├── [7.5G]  VH02001700_S6_R1_001.fastq.gz
            │   └── [17.8G]  VH02001700_S6_R2_001.fastq.gz
            ├── [19.0G]  SRR14730305
            │   └── [19.0G]  SRR14730305
            ├── SRR14730305
            │   ├── [8.4G]  VH02001702_S5_R1_001.fastq.gz
            │   └── [20.6G]  VH02001702_S5_R2_001.fastq.gz
            ├── [14.6G]  SRR14730306
            │   └── [14.6G]  SRR14730306
            ├── SRR14730306
            │   ├── [6.6G]  VH02001704_S4_R1_001.fastq.gz
            │   └── [16.0G]  VH02001704_S4_R2_001.fastq.gz
            ├── [21.5G]  SRR14730307
            │   └── [21.5G]  SRR14730307
            ├── SRR14730307
            │   ├── [9.6G]  VH02001708_S3_R1_001.fastq.gz
            │   └── [23.2G]  VH02001708_S3_R2_001.fastq.gz
            ├── [20.7G]  SRR14730308
            │   └── [20.7G]  SRR14730308
            ├── SRR14730308
            │   ├── [9.3G]  VH02001710_S2_R1_001.fastq.gz
            │   └── [22.1G]  VH02001710_S2_R2_001.fastq.gz
            ├── [15.8G]  SRR14730309
            │   └── [15.8G]  SRR14730309
            └── SRR14730309
                ├── [7.2G]  VH02001712_S1_R1_001.fastq.gz
                └── [16.9G]  VH02001712_S1_R2_001.fastq.gz

18 directories, 37 files
```


The `orig` directory contains the original fastq files. The fastq files are available for 10k and 100k subsamples in the `10k` and `100k` directories, respectively.

The `2-wells.fasta` file contains the barcodes for 2 wells.

## Test run

The pipeline can be run by creating a `params.yaml` file like this:

```yaml
param_list:
  - input_r1: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R1_001.fastq"
    input_r2: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R2_001.fastq"
    genomeDir: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.star.sparse"
    barcodesFasta: "gs://viash-hub-test-data/htrnaseq/v1/2-wells.fasta"
    id: sample_one
  - input_r1: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730302/VH02001614_S8_R1_001.fastq"
    input_r2: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730302/VH02001614_S8_R2_001.fastq"
    genomeDir: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.star.sparse"
    barcodesFasta: "gs://viash-hub-test-data/htrnaseq/v1/2-wells.fasta"
    id: sample_two
```

and then:

```bash
viash ns build --setup cb
nextflow run . -main-script target/nextflow/workflows/htrnaseq/main.nf \
  -profile docker \
  -c target/nextflow/workflows/htrnaseq/nextflow.config \
  -params-file params.yaml \
  -resume \
  --publish_dir output
```

Or, by running `src/workflows/htrnaseq/integration_test.sh`.
