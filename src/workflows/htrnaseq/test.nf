nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { htrnaseq } from targetDir + "/workflows/htrnaseq/main.nf"
include { check_eset } from targetDir + "/integration_test_components/htrnaseq/check_eset/main.nf"


params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  resources_test_file = file(params.resources_test)
  htrnaseq_ch = Channel.fromList([
      [
          id: "sample_one",
          input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
          input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
          genomeDir: resources_test_file.resolve("genomeDir/gencode.v41.star.sparse"),
          barcodesFasta: resources_test_file.resolve("360-wells-with-ids.fasta"),
          annotation: resources_test_file.resolve("genomeDir/gencode.v41.annotation.gtf.gz")
      ],
      [
          id: "sample_two",
          input_r1: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
          input_r2: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
          genomeDir: resources_test_file.resolve("genomeDir/gencode.v41.star.sparse"),
          barcodesFasta: resources_test_file.resolve("360-wells-with-ids.fasta"),
          annotation: resources_test_file.resolve("genomeDir/gencode.v41.annotation.gtf.gz")
      ]
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | htrnaseq.run(
        toState: {id, result, state -> result}
    )

  events_ch = htrnaseq_ch
    | toSortedList()
    | map {events ->
      assert events.size() == 2
      def ids = events.collect{it[0]}
      assert ids.toSet() == ["sample_one", "sample_two"].toSet()
      events.each {id, state ->
        assert state.fastq_output.each{it.isDirectory()}
        assert state.fastq_output.size() == 1
        assert state.star_output.each{it.isDirectory()}
        assert state.star_output.size() == 360
        assert state.f_data.isFile()
        assert state.p_data.isFile()
        assert state.html_report.isFile()
        assert state.run_params.isFile()
        assert state.star_qc_metrics.isFile()
        assert state.nrReadsNrGenesPerChrom.isFile()
        assert state.eset.isFile()
      }
      return events
    }

  check_eset_ch = htrnaseq_ch
    | check_eset.run(
        runIf: {id, state -> id == "sample_one"},
        args: [
          "expected_matrix": resources_test_file.resolve("expected_expressions.matrix")
        ],
        toState: [
            "eset": "eset",
            "star_output": "star_output"
        ]
    )
}


workflow test_wf2 {
  // Test the edge case where one of the barcodes has no reads
  resources_test_file = file(params.resources_test)
  input_ch = Channel.fromList([
      [
          id: "sample_one",
          input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
          input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
          genomeDir: resources_test_file.resolve("genomeDir/gencode.v41.star.sparse"),
          barcodesFasta: resources_test_file.resolve("2-wells-1-no-reads.fasta"),
          annotation: resources_test_file.resolve("genomeDir/gencode.v41.annotation.gtf.gz")
      ],
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | htrnaseq.run(
        toState: [
            "eset": "eset",
            "star_output": "star_output",
        ]
    )
}


workflow test_no_events {
  resources_test_file = file(params.resources_test)
  input_ch = Channel.fromList([])
    | map{ state -> [state.id, state] }
    | htrnaseq.run(
        toState: [
            "eset": "eset",
            "star_output": "star_output",
        ]
    )
    | toSortedList()
    | map {outputs ->
        assert outputs.size() == 0
    }
}


workflow test_concatenation {
  // Test concatenating samples
  // We give the two inputs the same sample ID, which causes them to be concatenated
  resources_test_file = file(params.resources_test)
  htrnaseq_ch = Channel.fromList([
      [
        id: "sample_one_run1",
        input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
        genomeDir: resources_test_file.resolve("genomeDir/gencode.v41.star.sparse"),
        barcodesFasta: resources_test_file.resolve("2-wells-1-no-reads.fasta"),
        annotation: resources_test_file.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        sample_id: "sample_one"
      ],
      [
        id: "sample_one_run2",
        input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
        genomeDir: resources_test_file.resolve("genomeDir/gencode.v41.star.sparse"),
        barcodesFasta: resources_test_file.resolve("2-wells-1-no-reads.fasta"),
        annotation: resources_test_file.resolve("genomeDir/gencode.v41.annotation.gtf.gz"),
        sample_id: "sample_one"
      ],
    ])
    | map{ state -> [state.id, state] }
    | view { "Input: $it" }
    | htrnaseq.run(
        toState: {id, result, state -> result}
    )

    events_ch = htrnaseq_ch
      | toSortedList()
      | map {events ->
        assert events.size() == 1
        def ids = events.collect{it[0]}
        assert ids.toSet() == ["sample_one"].toSet()
        events.each {id, state ->
            assert state.fastq_output.each{it.isDirectory()}
            assert state.fastq_output.size() == 2
            assert state.star_output.each{it.isDirectory()}
            assert state.star_output.size() == 2
            assert state.f_data.isFile()
            assert state.p_data.isFile()
            assert state.html_report.isFile()
            assert state.run_params.isFile()
            assert state.star_qc_metrics.isFile()
            assert state.nrReadsNrGenesPerChrom.isFile()
            assert state.eset.isFile()
        }
        return events
      }
}

