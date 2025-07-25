# htrnaseq v0.8.1

## Bug fixes

* Fix an issue where the FASTQ files from different samples on the same sequencing run would overwrite each other (PR #56).

## Under the hood

* Moved the test resources to their new location (PR #47).
* Add the package config (`_viash.yaml`) to every component's target dir. This makes introspection from, e.g. a `runner` workflow much more robust (PR #61)

## Minor changes

* Bump `biobox` and `craftbox` dependencies to versions `0.3.1` and `0.2.0`, respectively (PR #60).

# htrnaseq v0.8.0

## New functionality

* `save_params`: added a component to save workflow input parameters as yaml (PR #48).

* Added `run_params` parameter to `htrnaseq` and `runner` workflows in order to save the input parameters
  used for the workflow run (PR #48).
  
# htrnaseq v0.7.2

## Documentation

* Update README (PR #54)

# htrnaseq v0.7.1

## Bug fixes

* Bump viash version to `0.9.4`. This adds support for nextflow versions starting major version 25.01 and 
fixes an issue where an integer being passed to a argument with `type: double` resulted in an error (PR #51).

* `reporting`: updated default colour mapping (PR #50).

## Minor changes

* `create_report`: bump bioconductor version to 3.21 in order to accommodate R version 4.5 (PR #52).

# htrnaseq v0.7.0

## Breaking changes

The `runner` and `htrnaseq` workflow now output FASTQ files corresponding to the barcodes per input ID (per sequencing run). 
Previously, when multiple input folders or multiple input FASTQ files were provided
(for the `runner` and `htrnaseq` workflows respectively), the demultiplexed FASTQ files for these inputs were concatenated
and provided as output. For the `htrnaseq` workflow, reads can still be combined by using a newly added `sampleID` argument.
This means that two lists of FASTQ files can be provided for a single sample, and by assigning the same `sampleID`,
these reads will be joined. For example, with other arguments are left out for brevity:

```yaml
- id: sample1_run1
  input_r1: [sample_1_L001_1_R1.fastq, sample_1_L002_1_R1.fastq]
  input_r2: [sample_1_L001_1_R2.fastq, sample_1_L002_1_R2.fastq]
  sampleID: "sample_1"
- id: sample1_run2
  input_r1: [sample_1_L001_1_R1.fastq, sample_1_L002_1_R1.fastq]
  input_r2: [sample_1_L001_1_R2.fastq, sample_1_L002_1_R2.fastq]
  sampleID: "sample_1"
- id: sample_2
  input_r1: [sample_2_L001_1_R1.fastq, sample_2_L002_1_R1.fastq]
  input_r2: [sample_2_L001_1_R2.fastq, sample_2_L002_1_R2.fastq]
```

For the runner, concatenation of data across samples is automatically inferred. Previously, multiple IDs (events) could be
provided which were processed in parallel. This is no longer possible, as providing multiple will cause the matching
samples for these runs to be concatenated.


For example, the following old parameter yaml
```yaml
- id: run1
  input: ["run_folder_1/", run_folder_2/]
```
should now be provided as: 
```yaml
- id: run1
  input: "run_folder_1/"
- id: run2
  input: run_folder_2/
```

## Minor changes

* Updated viash to `0.9.2` (PR #49)

# htrnaseq v0.6.0

## Breaking changes

* `runner`: a subdirectory `data_processed` is now added to the output structure, in between
  the experiment ID and the directory with the workflow date and version (PR #45).

# htrnaseq v0.5.5

## New functionality

* Add `umi_length` parameter to the `runner` workflow (PR #46)

# htrnaseq v0.5.4

* Fix missing barcodes in the output from `generate_pool_statistics`, which caused an assertion error in `create_pdata`.
  In order to resolve the issue `generate_well_statistics` now outputs results for all chromosomes/scaffolds presented by
  the genome annotation, even when no reads were mapped to the regions in question. `generate_pool_statistics` will now
  remove regions from the output that have not at least one counts across all barcodes (PR #44).

# htrnaseq v0.5.3

## Bug fixes

* Fix `create_eset` component failing to create when one of the input samples has no counts (PR #43).

# htrnaseq v0.5.2

## Bug fixes

* `create_fdata`: remove duplicate entries from feature data (PR #41).

# htrnaseq v0.5.1

## Bug fixes

* `generate_well_statistics`: fix `ValueError` when an empty .bam file is provided as input (PR #40).
* `create_pdata`: avoid false positive `ValueError` for non-overlapping barcodes when input
  data contains empty (`NA`) values (PR #40).
  

# htrnaseq v0.5.0

## New functionality

* Added `ignore` parameter was added to the runner workflow in order to pass over certain input files
  from the input directory (PR #39).

# htrnaseq v0.4.0

## Breaking changes

An effort has been made to align the inputs for the `htrnaseq` and the mapping and demultiplexing of the wells, in order
simplify running these steps as seperate steps (PR #37).
  * Changes to the `parallel_map` component:
    - The `barcode` argument has been renamed to `barcodesFasta` and the provided 
      value for this argument must now be single FASTA file instead of a list of barcodes.
    - The filenames for the provided FASTQ files must now conform to the format `{name}_R(1|2).fasta`,
      where `{name}` is the well identifiers. The well identifiers correspond to the headers
      of the FASTA file containing the barcodes (up untill the first whitespace).
      Forward and reverse FASTQ files must still be provided in pairs, meaning that the order of
      files provided to `input_r1` and `input_r2` remains important.
    - The requirement for equal number of barcodes and FASTQ pairs to be provided has been dropped.
      Instead, the barcodes provided with `barcodesFasta` are matched to the input FASTQ files by comparing
      the header of the FASTA records to the file names of the provided FASTQ input files. Each barcode must
      match exactly one FASTQ input pair (forward and reverse reads), but FASTQ files that were not matched to any
      barcode are not processed. Basically, the barcodes fasta can now act as a filter for the FASTQ files to be mapped.
  * The `utils/groupWells` workflow has been removed.
  * `parallel_map_wf` has been removed as its functionality is now incomporated into the `parallel_map` component. 
  * The `pool`, `well_id`, `barcode`, `lane`, `pair_end` and `n_wells` output arguments have been dropped from the 
     `well_demultiplexing` workflow. This workflow now only outputs a list of demultiplexed FASTQ files.
  * A `well_metadata` workflow has been implemented that extracts the metadata that is no longer output by the `well_demultiplexing`
    workflow from the demultiplexed files and the barcodes FASTA.

## New functionality

* Multiple input directories can not be provided. The input reads from these from these directories
  will be joined per barcode before mapping. This is useful when data has been generated using
  multiple sequencing runs in order to increase sequencing depth (PR #38).

# htrnaseq v0.3.0

## New functionality

* Added `umi_length` argument (PR #27).
* Added `runner` workflow (PR #26, see below)

## `runner` workflow

* Removed `wellBarcodesLength` from `parallel_map` workflow (PR #27).

## Major changes

A runner workflows has been added, providing two additional features:

1. Start from an input directory containing fastq files rather than a list of input fastq pairs.
2. Improve the output of the workflow

### Input directory

It is now possible to specify a single `--input <basedir>` directory as input and the runner will extract the fastq file pairs. An error will be raised if the filename processing leads to errors.

### Output

The runner provides a complete different approach to output. A couple of things are important here:

- Output is split up in 2 parts:

    1. The well-demultiplexed fastq files (`--fastq_publish_dir`)
    2. All the other results of the workflow (`--results_publish_dir`)

- The well-demultiplexed fastq file are stored under `--fastq_publish_dir` according to the following format:

    ```
    $fastq_publish_dir/$id/<date-time>_htrnaseq_<version>/$plate_$lane/<well_id>_R1/2_001.fastq
    ```

- The other results are stored under `--results_publish_dir` according to the following format:

    ```
    $results_publish_dir/$project_id/$experiment_id/<date-time>_htrnaseq_<version>/
    ```

    This is an example listing of this directory:

    ```
    esets
    fData
    nrReadsNrGenesPerChrom
    pData
    report.html
    star_output
    starLogs
    ```

This output structure can be circumvented by using the `--output_dir` option, which will store all output in a single directory.

1. Using the `htrnaseq` workflow directory rather than the `runner` interface
2. Using the argument `--plain_output` with the `runner`. fastq files and other results will still be published in their respective directories, but not in a directory hierarchy as described above.

## Minor changes

* Use `v0.2.0` version of cutadapt instead of `main` (PR #23).
* Use `v0.3.0` version of cutadapt
* Bump viash to 0.9.1 (PR #31).
* `create_eset`: Update base container image, `R` version and all dependencies
  to newer versions (PR #28).

# htrnaseq v0.2.0

# New functionality

* Make sure that the Well ID matches the required format (PR #22 and PR #21). 

# htrnaseq v0.1.0

Initial release
