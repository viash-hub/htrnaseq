include { well_demultiplex } from params.rootDir + "/target/nextflow/workflows/well_demultiplex/main.nf"
include { check_cutadapt_output } from params.rootDir + "/target/nextflow/integration_test_components/well_demultiplexing/check_cutadapt_output/main.nf"

params.resources_test =  "gs://viash-hub-test-data/htrnaseq/v1/"

workflow test_wf {
  resources_test_file = file(params.resources_test)
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
        barcodesFasta: resources_test_file.resolve("2-wells-with-ids.fasta"),
      ],
      [
        id: "SRR14730302",
        input_r1: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
        input_r2: resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
        barcodesFasta:  resources_test_file.resolve("2-wells-with-ids.fasta"),
      ],
    ])
    | map { state -> [ state.id, state ] }
    | well_demultiplex.run(
      fromState: { id, state ->
        [
          input_r1: state.input_r1,
          input_r2: state.input_r2,
          barcodesFasta: state.barcodesFasta,
        ]
      },
      toState: { id, output, state ->
        output }
    )
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 4 : "2 samples, each with 2 barcodes"
    }
    | map {output ->
      def ids = output.collect{it[0]}
      def states = output.collect{it[1]}
      def new_state = [
        "ids": ids,
        "fastq_r1": states.collect{it.output_r1}.flatten(),
        "fastq_r2": states.collect{it.output_r2}.flatten()
      ]
      ["integration_test_check", new_state]
    }
    | check_cutadapt_output.run(
      fromState: {id, state -> state}
    )
}


workflow test_wf2 {
  resources_test_file = file(params.resources_test)
  output_ch = Channel.fromList([
      [
        id: "SRR14730301",
        input_r1:
          [
            resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R1_001.fastq"),
            resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R1_001.fastq"),
          ],
        input_r2:
          [
            resources_test_file.resolve("100k/SRR14730301/VH02001612_S9_R2_001.fastq"),
            resources_test_file.resolve("100k/SRR14730302/VH02001614_S8_R2_001.fastq"),
          ],
        barcodesFasta: resources_test_file.resolve("2-wells-with-ids.fasta"),
      ],
    ])
    | map { state -> [ state.id, state ] }
    | well_demultiplex.run(
      fromState: { id, state ->
        [
          input_r1: state.input_r1,
          input_r2: state.input_r2,
          barcodesFasta: state.barcodesFasta,
        ]
      },
      toState: { id, output, state ->
        output }
    )
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      "Output: $output"
    }
    | toSortedList()
    | view { output ->
      assert output.size() == 2 : "1 samples, and two barcodes"
    }
}
