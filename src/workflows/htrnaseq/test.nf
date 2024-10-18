nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { htrnaseq } from targetDir + "/workflows/htrnaseq/main.nf"
include { check_eset } from targetDir + "/integration_test_components/htrnaseq/check_eset/main.nf"


params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  input_ch = Channel.fromList([
      [
          id: "sample_one",
          input_r1: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R1_001.fastq",
          input_r2: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730301/VH02001612_S9_R2_001.fastq",
          genomeDir: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.star.sparse",
          barcodesFasta: "gs://viash-hub-test-data/htrnaseq/v1/360-wells.fasta",
          annotation: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.annotation.gtf.gz"
      ],
      [
          id: "sample_two",
          input_r1: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730302/VH02001614_S8_R1_001.fastq",
          input_r2: "gs://viash-hub-test-data/htrnaseq/v1/100k/SRR14730302/VH02001614_S8_R2_001.fastq",
          genomeDir: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.star.sparse",
          barcodesFasta: "gs://viash-hub-test-data/htrnaseq/v1/360-wells.fasta",
          annotation: "gs://viash-hub-test-data/htrnaseq/v1/genomeDir/gencode.v41.annotation.gtf.gz"
      ]
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | htrnaseq.run(
        toState: [
            "eset": "eset",
            "star_output": "star_output",
        ]
    )
    | check_eset.run(
        runIf: {id, state -> id == "sample_one"},
        toState: [
            "eset": "eset",
            "star_output": "star_output"
        ]
    )
}