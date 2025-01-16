# demultiplex v0.x.x

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

# demultiplex v0.2.0

# New functionality

* Make sure that the Well ID matches the required format (PR #22 and PR #21). 

# demultiplex v0.1.0

Initial release
